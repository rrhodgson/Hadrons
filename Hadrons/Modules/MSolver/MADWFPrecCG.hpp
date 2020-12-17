/*
 * MADWFPrecCG.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MSolver_MADWFPrecCG_hpp_
#define Hadrons_MSolver_MADWFPrecCG_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MSolver/Guesser.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Mixed precision schur red-black preconditioned CG             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class MADWFPrecCGPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MADWFPrecCGPar,
                                    std::string , innerAction,
                                    std::string , outerAction,
                                    unsigned int, maxInnerIteration,
                                    unsigned int, maxOuterIteration,
                                    double      , residual,
                                    std::string , eigenPack,
                                    std::string , gaugeOuter,
                                    std::string , gaugeInner);
};

template <typename FImplInner, typename FImplOuter, int nBasis, typename GImpl>
class TMADWFPrecCG: public Module<MADWFPrecCGPar>
{
public:
    FERM_TYPE_ALIASES(FImplInner, Inner);
    FERM_TYPE_ALIASES(FImplOuter, Outer);
    SOLVER_TYPE_ALIASES(FImplOuter,);
    // GAUGE_TYPE_ALIASES(GImpl,);
    typedef HADRONS_DEFAULT_SCHUR_OP<FMatInner, FermionFieldInner> SchurFMatInner;
    typedef HADRONS_DEFAULT_SCHUR_OP<FMatOuter, FermionFieldOuter> SchurFMatOuter;
private:
    template <typename Field>
    class OperatorFunctionWrapper: public OperatorFunction<Field>
    {
    public:
        using OperatorFunction<Field>::operator();
        OperatorFunctionWrapper(LinearFunction<Field> &fn): fn_(fn) {};
        virtual ~OperatorFunctionWrapper(void) = default;
        virtual void operator()(LinearOperatorBase<Field> &op, 
                                const Field &in, Field &out)
        {
            fn_(in, out);
        }
    private:
        LinearFunction<Field> &fn_;
    };
public:
    // constructor
    TMADWFPrecCG(const std::string name);
    // destructor
    virtual ~TMADWFPrecCG(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

// MODULE_REGISTER_TMP(MADWFCG, ARG(TMADWFPrecCG<ZFIMPLD, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS, GIMPL>), MSolver);
MODULE_REGISTER_TMP(MixedPrecMADWFCG, ARG(TMADWFPrecCG<ZFIMPLD, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS, GIMPL>), MSolver);

/******************************************************************************
 *                 TMADWFPrecCG implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis, typename GImpl>
TMADWFPrecCG<FImplInner, FImplOuter, nBasis, GImpl>
::TMADWFPrecCG(const std::string name)
: Module<MADWFPrecCGPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis, typename GImpl>
std::vector<std::string> TMADWFPrecCG<FImplInner, FImplOuter, nBasis, GImpl>
::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImplInner, typename FImplOuter, int nBasis, typename GImpl>
std::vector<std::string> TMADWFPrecCG<FImplInner, FImplOuter, nBasis, GImpl>
::getReference(void)
{
    std::vector<std::string> ref = {par().innerAction, par().outerAction, par().gaugeInner, par().gaugeOuter};
    
    if (!par().eigenPack.empty())
    {
        ref.push_back(par().eigenPack);
    }
    
    return ref;
}

template <typename FImplInner, typename FImplOuter, int nBasis, typename GImpl>
std::vector<std::string> TMADWFPrecCG<FImplInner, FImplOuter, nBasis, GImpl>
::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_subtract"};

    return out;
}

struct CGincreaseTol : public MADWFinnerIterCallbackBase{
  ConjugateGradient<LatticeFermionD> &cg_inner;  
  RealD outer_resid;

  CGincreaseTol(ConjugateGradient<LatticeFermionD> &cg_inner,
        RealD outer_resid): cg_inner(cg_inner), outer_resid(outer_resid){}
  
  void operator()(const RealD current_resid){
    std::cout << "CGincreaseTol with current residual " << current_resid << " changing inner tolerance " << cg_inner.Tolerance << " -> ";
    while(cg_inner.Tolerance < current_resid) cg_inner.Tolerance *= 2;    
    //cg_inner.Tolerance = outer_resid/current_resid;
    std::cout << cg_inner.Tolerance << std::endl;
  }
};


// setup ///////////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis, typename GImpl>
void TMADWFPrecCG<FImplInner, FImplOuter, nBasis, GImpl>
::setup(void)
{
    LOG(Message) << "Setting up Schur red-black preconditioned MADWF "
                 << "CG for inner/outer action '" << par().innerAction 
                 << "'/'" << par().outerAction << "', residual "
                 << par().residual << ", and maximum inner/outer iteration " 
                 << par().maxInnerIteration << "/" << par().maxOuterIteration
                 << std::endl;


    auto Ls_outer        = env().getObjectLs(par().outerAction);
    auto Ls_inner        = env().getObjectLs(par().innerAction);

    auto &Umu_outer = envGet(GaugeFieldOuter, par().gaugeOuter);
    auto &Umu_inner = envGet(GaugeFieldInner, par().gaugeInner);

    auto &omat     = envGet(FMatOuter, par().outerAction);
    // auto &imat     = envGet(FMatInner, par().innerAction);

    MobiusFermionD &D_outer  = envGetDerived(FMatOuter,MobiusFermionD,par().outerAction);
    ZMobiusFermionD &D_inner = envGetDerived(FMatInner,ZMobiusFermionD,par().innerAction);

    double residual = par().residual;


std::cout << "Setup Mob action" << std::endl;

    auto &g4_outer   = *envGetGrid(FImplOuter::FermionField);
    auto &grb4_outer = *envGetRbGrid(FImplOuter::FermionField);

    auto &g5_outer   = *envGetGrid(FImplOuter::FermionField, Ls_outer);
    auto &grb5_outer = *envGetRbGrid(FImplOuter::FermionField, Ls_outer);
    MobiusFermionD D_outer_loc(Umu_outer, g5_outer, grb5_outer, g4_outer, grb4_outer, D_outer.mass, 1.8, D_outer._b, D_outer._c);


std::cout << "Setup zMob action" << std::endl;

    auto &g4_inner   = *envGetGrid(FImplInner::FermionField);
    auto &grb4_inner = *envGetRbGrid(FImplInner::FermionField);

    auto &g5_inner   = *envGetGrid(FImplInner::FermionField, Ls_inner);
    auto &grb5_inner = *envGetRbGrid(FImplInner::FermionField, Ls_inner);
    std::vector<ComplexD> gamma(D_inner._gamma.size());
    for (int i=0; i<D_inner._gamma.size(); i++)
        gamma[i] = D_inner._gamma[i];
    ZMobiusFermionD D_inner_loc(Umu_inner, g5_inner, grb5_inner, g4_inner, grb4_inner, D_inner.mass, 1.8, gamma, D_inner._b, D_inner._c);


    auto makeSolver = [&D_outer_loc, &D_inner_loc, &omat, Ls_outer, Ls_inner, &Umu_outer, &Umu_inner, residual, this] (bool subGuess) mutable
    {
        return [&D_outer_loc, &D_inner_loc, &omat, Ls_outer, Ls_inner, &Umu_outer, &Umu_inner, subGuess, residual, this]
        (FermionFieldOuter &sol, const FermionFieldOuter &source) mutable
        {


double resid_outer = residual;
double resid_inner = residual;


std::cout << "Setup source" << std::endl;


  // std::vector<int> seeds4({1, 2, 3, 4});

  // GridParallelRNG RNG4(Umu.Grid());
  // RNG4.SeedFixedIntegers(seeds4);

  LatticeFermionD src4(Umu_outer.Grid());
  //random(RNG4,src4);
  D_outer_loc.ExportPhysicalFermionSource(source,src4);

  
  //Solve using a regular even-odd preconditioned CG for the Hermitian operator
  //M y = x
  //Mprec y'_o = x'_o     where Mprec = Doo - Doe Dee^-1 Deo    and  x'_o = -Doe Dee^-1 x_e + x_o
  //y_o = y'_o

  //(Mprec^dag Mprec) y'_o = Mprec^dag x'_o 
  //y'_o = (Mprec^dag Mprec)^-1 Mprec^dag x'_o 


std::cout << "Setup Solvers" << std::endl;

  //We can get Mprec^dag x'_o from x_o  from SchurRedBlackDiagMooeeSolve::RedBlackSource
  ConjugateGradient<LatticeFermionD> CG_outer(resid_outer, par().maxInnerIteration);

  GridStopWatch CGTimer;


  typedef PauliVillarsSolverFourierAccel<LatticeFermionD, LatticeGaugeFieldD> PVtype;
  PVtype PV_outer(Umu_outer, CG_outer);

  ConjugateGradient<LatticeFermionD> CG_inner(resid_inner, par().maxInnerIteration, 0);

  CGincreaseTol update(CG_inner, resid_outer);

  SchurRedBlackDiagTwoSolve<LatticeFermionD> SchurSolver_inner(CG_inner);


  LatticeFermionD result_MADWF(omat.FermionGrid());
  result_MADWF = Zero();


std::cout << "Run MADWF" << std::endl;

  if (par().eigenPack.empty()) {
    ZeroGuesser<LatticeFermionD> Guess;
    MADWF<MobiusFermionD, ZMobiusFermionD, PVtype, SchurRedBlackDiagTwoSolve<LatticeFermionD>, ZeroGuesser<LatticeFermionD> > 
                  madwf(D_outer_loc, D_inner_loc, PV_outer, SchurSolver_inner, Guess, resid_outer, par().maxOuterIteration, &update);

    CGTimer.Start();
    madwf(src4, result_MADWF);
    CGTimer.Stop();
  } else {
    auto guesserPt = makeGuesser<FImplInner, nBasis>(par().eigenPack);
    MADWF<MobiusFermionD, ZMobiusFermionD, PVtype, SchurRedBlackDiagTwoSolve<LatticeFermionD>, LinearFunction<LatticeFermionD> >  
                  madwf(D_outer_loc, D_inner_loc, PV_outer, SchurSolver_inner, *guesserPt, resid_outer, par().maxOuterIteration, &update);

    CGTimer.Start();
    madwf(src4, result_MADWF);
    CGTimer.Stop();
  }
  

  std::cout << GridLogMessage << "Total MADWF time : " << CGTimer.Elapsed()
            << std::endl;


  std::cout << GridLogMessage << "######## Dhop calls summary" << std::endl;
  D_outer_loc.Report();

    sol = result_MADWF;



        };
    };
    auto solver = makeSolver(false);
    envCreate(Solver, getName(), Ls_outer, solver, omat);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls_outer, solver_subtract, omat);
}






// execution ///////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis, typename GImpl>
void TMADWFPrecCG<FImplInner, FImplOuter, nBasis, GImpl>
::execute(void)
{
/*
std::string config_file = par().gaugefile; //"/home/dp008/dp008/dc-hodg1/Gauge_Confs/16^3/ckpoint_lat.IEEE64BIG.1100";
bool load_config = true;
if (config_file.empty())
    load_config = false;

double mass = 0.01;

double b_plus_c_inner = 1.0;
int Ls_inner = 24;

double b_plus_c_outer = 2.0;
int Ls_outer = 10;

double lambda_max = 1.42;

double resid_outer = 1e-8;
double resid_inner = 1e-8;

    RealD bmc = 1.0; //use Shamir kernel
  std::vector<ComplexD> gamma_inner;

//   std::cout << "Compute parameters" << std::endl;
// Approx::computeZmobiusGamma(gamma_inner, b_plus_c_inner, Ls_inner, b_plus_c_outer, Ls_outer, lambda_max);


  // std::cout << "Accept parameters" << std::endl;

    // For 24^3 w/ deflation
    // gamma_inner.push_back(ComplexD(1.38355,-0));
    // gamma_inner.push_back(ComplexD(1.15775,-0));
    // gamma_inner.push_back(ComplexD(0.791841,-0));
    // gamma_inner.push_back(ComplexD(0.491647,-0));
    // gamma_inner.push_back(ComplexD(0.289819,-0));
    // gamma_inner.push_back(ComplexD(0.1744,0.0227481));
    // gamma_inner.push_back(ComplexD(0.145524,0.0736917));
    // gamma_inner.push_back(ComplexD(0.0984094,0.13761));
    // gamma_inner.push_back(ComplexD(0.0984094,-0.13761));
    // gamma_inner.push_back(ComplexD(0.145524,-0.0736917));
    // gamma_inner.push_back(ComplexD(0.1744,-0.0227481));
    // gamma_inner.push_back(ComplexD(0.218611,-0));
    // gamma_inner.push_back(ComplexD(0.379624,-0));
    // gamma_inner.push_back(ComplexD(0.628981,-0));
    // gamma_inner.push_back(ComplexD(0.97407,-0));
    // gamma_inner.push_back(ComplexD(1.33346,-0));



    // // For 16^3 or 48^3 w/ no deflation
    gamma_inner.push_back(ComplexD(1.458064389850479e+00,-0.000000000000000e+00));
    gamma_inner.push_back(ComplexD(1.182313183893475e+00,-0.000000000000000e+00));
    gamma_inner.push_back(ComplexD(8.309511666859551e-01,-0.000000000000000e+00));
    gamma_inner.push_back(ComplexD(5.423524091567911e-01,-0.000000000000000e+00));
    gamma_inner.push_back(ComplexD(3.419850204537295e-01,-0.000000000000000e+00));
    gamma_inner.push_back(ComplexD(2.113790261902896e-01,-0.000000000000000e+00));
    gamma_inner.push_back(ComplexD(1.260742995029118e-01,-0.000000000000000e+00));
    gamma_inner.push_back(ComplexD(9.901366519626265e-02,-0.000000000000000e+00));
    gamma_inner.push_back(ComplexD(6.863249884465925e-02,5.506585308274019e-02));
    gamma_inner.push_back(ComplexD(6.863249884465925e-02,-5.506585308274019e-02));

  
  std::cout << "gamma:\n";
  for(int s=0;s<Ls_inner;s++) std::cout << s << " " << gamma_inner[s] << std::endl;


  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()),
      GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);


  GridCartesian* FGrid_outer = SpaceTimeGrid::makeFiveDimGrid(Ls_outer, UGrid);
  GridCartesian* FGrid_inner = SpaceTimeGrid::makeFiveDimGrid(Ls_inner, UGrid);

  GridRedBlackCartesian* FrbGrid_outer = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls_outer, UGrid);
  GridRedBlackCartesian* FrbGrid_inner = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls_inner, UGrid);


  std::vector<int> seeds4({1, 2, 3, 4});
  std::vector<int> seeds5({5, 6, 7, 8});

  GridParallelRNG RNG5_outer(FGrid_outer);
  RNG5_outer.SeedFixedIntegers(seeds5);

  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);


  LatticeFermionD src4(UGrid); random(RNG4,src4);


  LatticeFermionD result_outer(FGrid_outer);
  result_outer = Zero();
  LatticeGaugeFieldD Umu(UGrid);

  if(load_config){
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, config_file);

    for(int i=0;i<Nd;i++){
      assert(header.dimension[i] == GridDefaultLatt()[i]);
    }
  }else{    
    // SU<Nc>::HotConfiguration(RNG4, Umu); // Random Gauge
    SU<Nc>::ColdConfiguration(RNG4, Umu); // Unit Gauge
  }
    
  std::cout << GridLogMessage << "Lattice dimensions: " << GridDefaultLatt()
            << "   Ls: " << Ls_outer << std::endl;

  RealD M5 = 1.8;

  RealD b_outer = (b_plus_c_outer + bmc)/2.;
  RealD c_outer = (b_plus_c_outer - bmc)/2.;

  RealD b_inner = (b_plus_c_inner + bmc)/2.;
  RealD c_inner = (b_plus_c_inner - bmc)/2.;

  MobiusFermionD D_outer(Umu, *FGrid_outer, *FrbGrid_outer, *UGrid, *UrbGrid, mass, M5, b_outer, c_outer);
  ZMobiusFermionD D_inner(Umu, *FGrid_inner, *FrbGrid_inner, *UGrid, *UrbGrid, mass, M5, gamma_inner, b_inner, c_inner);


  // LatticeFermionD src_outer(FGrid_outer);
  // D_outer.ImportPhysicalFermionSource(src4,src_outer); //applies D_- 

  //Solve using a regular even-odd preconditioned CG for the Hermitian operator
  //M y = x
  //Mprec y'_o = x'_o     where Mprec = Doo - Doe Dee^-1 Deo    and  x'_o = -Doe Dee^-1 x_e + x_o
  //y_o = y'_o

  //(Mprec^dag Mprec) y'_o = Mprec^dag x'_o 
  //y'_o = (Mprec^dag Mprec)^-1 Mprec^dag x'_o 

  //We can get Mprec^dag x'_o from x_o  from SchurRedBlackDiagMooeeSolve::RedBlackSource
  ConjugateGradient<LatticeFermionD> CG_outer(resid_outer, 10000);

  GridStopWatch CGTimer;
  
  typedef PauliVillarsSolverFourierAccel<LatticeFermionD, LatticeGaugeFieldD> PVtype;
  PVtype PV_outer(Umu, CG_outer);

  ConjugateGradient<LatticeFermionD> CG_inner(resid_inner, 10000, 0);

  CGincreaseTol update(CG_inner, resid_outer);

  SchurRedBlackDiagTwoSolve<LatticeFermionD> SchurSolver_inner(CG_inner);

  ZeroGuesser<LatticeFermion> Guess;
  MADWF<MobiusFermionD, ZMobiusFermionD, PVtype, SchurRedBlackDiagTwoSolve<LatticeFermionD>, ZeroGuesser<LatticeFermion> > madwf(D_outer, D_inner, PV_outer, SchurSolver_inner, Guess, resid_outer, 100, &update);
  
  LatticeFermionD result_MADWF(FGrid_outer);
  result_MADWF = Zero();

  CGTimer.Start();
  madwf(src4, result_MADWF);
  CGTimer.Stop();

  LatticeFermionD result_o_MADWF(FrbGrid_outer);
  pickCheckerboard(Odd, result_o_MADWF, result_MADWF);

  std::cout << GridLogMessage << "Total MADWF time : " << CGTimer.Elapsed()
            << std::endl;

  std::cout << GridLogMessage << "######## Dhop calls summary" << std::endl;
  D_outer.Report();
*/
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_MADWFPrecCG_hpp_
