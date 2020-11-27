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
                                    std::string , gauge,
                                    std::string , gaugefile);
};

template <typename FImplInner, typename FImplOuter, int nBasis, typename GImpl>
class TMADWFPrecCG: public Module<MADWFPrecCGPar>
{
public:
    FERM_TYPE_ALIASES(FImplInner, Inner);
    FERM_TYPE_ALIASES(FImplOuter, Outer);
    SOLVER_TYPE_ALIASES(FImplOuter,);
    GAUGE_TYPE_ALIASES(GImpl,);
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

// MODULE_REGISTER_TMP(MADWFPrecCG, 
    // ARG(TMADWFPrecCG<ZFIMPLF, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
MODULE_REGISTER_TMP(ZMADWFPrecCG, ARG(TMADWFPrecCG<ZFIMPLD, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS, GIMPL>), MSolver);

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
    std::vector<std::string> ref = {par().innerAction, par().outerAction, par().gauge};
    
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

    auto &Umu = envGet(GaugeField, par().gauge);

    auto &omat     = envGet(FMatOuter, par().outerAction);
    // auto &imat     = envGet(FMatInner, par().innerAction);

    MobiusFermionD &D_outer  = envGetDerived(FMatOuter,MobiusFermionD,par().outerAction);
    ZMobiusFermionD &D_inner = envGetDerived(FMatInner,ZMobiusFermionD,par().innerAction);

    // auto &g4   = *envGetGrid(FImplOuter::FermionField);
    // auto &grb4 = *envGetRbGrid(FImplOuter::FermionField);
    // auto &g5   = *envGetGrid(FImplOuter::FermionField, Ls_outer);
    // auto &grb5 = *envGetRbGrid(FImplOuter::FermionField, Ls_outer);



    auto guesserPt = makeGuesser<FImplInner, nBasis>(par().eigenPack); 

    double residual = par().residual;

    auto makeSolver = [&D_outer, &D_inner, &omat, &Umu, guesserPt, residual, this] (bool subGuess) mutable
    {
        return [&D_outer, &D_inner, &omat, &Umu, guesserPt, subGuess, residual, this]
        (FermionFieldOuter &sol, const FermionFieldOuter &source) mutable
        {


    auto &g4   = *envGetGrid(FImplOuter::FermionField);
    auto &grb4 = *envGetRbGrid(FImplOuter::FermionField);

    auto &g5_outer   = *envGetGrid(FImplOuter::FermionField, 24);
    auto &grb5_outer = *envGetRbGrid(FImplOuter::FermionField, 24);
    MobiusFermionD D_outer_loc(Umu, g5_outer, grb5_outer, g4, grb4, D_outer.mass, 1.8, D_outer._b, D_outer._c);

    auto &g5_inner   = *envGetGrid(FImplInner::FermionField, 10);
    auto &grb5_inner = *envGetRbGrid(FImplInner::FermionField, 10);
    std::vector<ComplexD> gamma(D_inner._gamma.size());
    for (int i=0; i<D_inner._gamma.size(); i++)
        gamma[i] = D_inner._gamma[i];
    ZMobiusFermionD D_inner_loc(Umu, g5_inner, grb5_inner, g4, grb4, D_inner.mass, 1.8, gamma, D_inner._b, D_inner._c);

std::cout << "mass = " << D_inner_loc.mass << std::endl;
  std::cout << "b = " << D_inner_loc._b << std::endl;
  std::cout << "c = " << D_inner_loc._c << std::endl;
  std::cout << "_zolo_hi = " << D_inner_loc._zolo_hi << std::endl;
  std::cout << "gamma = " << std::endl;
  for (int i=0; i<D_inner_loc._gamma.size(); i++) {
    std::cout << "\t" << D_inner_loc._gamma[i] << std::endl;
  }
  std::cout << "omega = " << std::endl;
  for (int i=0; i<D_inner_loc.omega.size(); i++) {
    std::cout << "\t" << D_inner_loc.omega[i] << std::endl;
  }
    std::cout << "bs = " << std::endl;
  for (int i=0; i<D_inner_loc.bs.size(); i++) {
    std::cout << "\t" << D_inner_loc.bs[i] << std::endl;
  }
  std::cout << "cs = " << std::endl;
  for (int i=0; i<D_inner_loc.cs.size(); i++) {
    std::cout << "\t" << D_inner_loc.cs[i] << std::endl;
  }
  std::cout << "as = " << std::endl;
  for (int i=0; i<D_inner_loc.as.size(); i++) {
    std::cout << "\t" << D_inner_loc.as[i] << std::endl;
  }


double resid_outer = residual;
double resid_inner = residual;


  std::vector<int> seeds4({1, 2, 3, 4});

  GridParallelRNG RNG4(Umu.Grid());
  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermionD src4(Umu.Grid()); random(RNG4,src4);

// SU<Nc>::ColdConfiguration(RNG4, Umu); // Unit Gauge
  
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
  // MADWF<MobiusFermionD, ZMobiusFermionD, PVtype, SchurRedBlackDiagTwoSolve<LatticeFermionD>, ZeroGuesser<LatticeFermion> > 
  //               madwf(D_outer_loc, D_inner_loc, PV_outer, SchurSolver_inner, Guess, resid_outer, 10, &update);
  MADWF<MobiusFermionD, ZMobiusFermionD, PVtype, SchurRedBlackDiagTwoSolve<LatticeFermionD>, LinearFunction<LatticeFermion> >  
                madwf(D_outer_loc, D_inner_loc, PV_outer, SchurSolver_inner, *guesserPt, resid_outer, 10, &update);
  

  LatticeFermionD result_MADWF(omat.FermionGrid());
  result_MADWF = Zero();

  CGTimer.Start();
  madwf(src4, result_MADWF);
  CGTimer.Stop();

  LatticeFermionD result_o_MADWF(omat.FermionRedBlackGrid());
  pickCheckerboard(Odd, result_o_MADWF, result_MADWF);

  std::cout << GridLogMessage << "Total MADWF time : " << CGTimer.Elapsed()
            << std::endl;



  std::cout << GridLogMessage << "######## Dhop calls summary" << std::endl;
  D_outer_loc.Report();







































            // std::cout << "Setup lattices" << std::endl;

            // LatticeGaugeFieldD Umu(U);

            // LatticeFermionD src_outer(source.Grid());

            // LatticeFermionD src4(U.Grid()); // TODO: Get this from source
            // GridParallelRNG RNG4(U.Grid());
            // std::vector<int> seeds4({1, 2, 3, 4});
            // RNG4.SeedFixedIntegers(seeds4);
            // random(RNG4,src4)


            // // D_outer.ImportPhysicalFermionSource(src4,src_outer); //applies D_-
            // // D_outer.ImportUnphysicalFermion(source,src4);
            // // src4 = source;

            // LatticeFermionD result_MADWF(source.Grid());
            // result_MADWF = Zero();


            // std::cout << "Setup CG and solvers" << std::endl;

            // ConjugateGradient<LatticeFermionD> CG_outer(residual, 10000);
            // ConjugateGradient<LatticeFermionD> CG_inner(residual, 10000, 0);
            // CGincreaseTol update(CG_inner, residual);

            // typedef PauliVillarsSolverFourierAccel<LatticeFermionD, LatticeGaugeFieldD> PVtype;
            // PVtype PV_outer(Umu, CG_outer);

            // typedef SchurRedBlackDiagTwoSolve<LatticeFermionD> SchurSolverType;
            // SchurSolverType SchurSolver_inner(CG_inner);


            // std::cout << "Setup MADWF" << std::endl;

            // MADWF<MobiusFermionD, ZMobiusFermionD, PVtype, SchurSolverType, LinearFunction<LatticeFermion> > 
            //         madwf(D_outer, D_inner, PV_outer, SchurSolver_inner, *guesserPt, residual, 100, &update);


            // std::cout << "Run MADWF" << std::endl;

            // madwf(src4, result_MADWF);

            // sol = result_MADWF;


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

std::string config_file = par().gaugefile; //"/home/dp008/dp008/dc-hodg1/Gauge_Confs/16^3/ckpoint_lat.IEEE64BIG.1100";
bool load_config = true;
if (config_file.empty())
    load_config = false;

double mass = 0.01;

double b_plus_c_inner = 1.0;
int Ls_inner = 10;

double b_plus_c_outer = 2.0;
int Ls_outer = 24;

double lambda_max = 1.42;

double resid_outer = 1e-8;
double resid_inner = 1e-8;

    RealD bmc = 1.0; //use Shamir kernel
  std::vector<ComplexD> gamma_inner;

//   std::cout << "Compute parameters" << std::endl;
// Approx::computeZmobiusGamma(gamma_inner, b_plus_c_inner, Ls_inner, b_plus_c_outer, Ls_outer, lambda_max);


  std::cout << "Accept parameters" << std::endl;

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


  // std::cout << src4 << std::endl;

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

  std::cout << "mass = " << D_inner.mass << std::endl;
  std::cout << "b = " << D_inner._b << std::endl;
  std::cout << "c = " << D_inner._c << std::endl;
  std::cout << "_zolo_hi = " << D_inner._zolo_hi << std::endl;
  std::cout << "gamma = " << std::endl;
  for (int i=0; i<D_inner._gamma.size(); i++) {
    std::cout << "\t" << D_inner._gamma[i] << std::endl;
  }
  std::cout << "omega = " << std::endl;
  for (int i=0; i<D_inner.omega.size(); i++) {
    std::cout << "\t" << D_inner.omega[i] << std::endl;
  }
  std::cout << "bs = " << std::endl;
  for (int i=0; i<D_inner.bs.size(); i++) {
    std::cout << "\t" << D_inner.bs[i] << std::endl;
  }
  std::cout << "cs = " << std::endl;
  for (int i=0; i<D_inner.cs.size(); i++) {
    std::cout << "\t" << D_inner.cs[i] << std::endl;
  }
  std::cout << "as = " << std::endl;
  for (int i=0; i<D_inner.as.size(); i++) {
    std::cout << "\t" << D_inner.as[i] << std::endl;
  }

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
  // SchurRedBlackDiagTwoSolve<LatticeFermionD> SchurSolver_outer(CG_outer);
  
  // LatticeFermionD tmp_e_outer(FrbGrid_outer);
  // LatticeFermionD src_o_outer(FrbGrid_outer);
  // SchurSolver_outer.RedBlackSource(D_outer, src_outer, tmp_e_outer, src_o_outer);
  
  // LatticeFermionD result_o_outer(FrbGrid_outer);
  // result_o_outer = Zero();

  GridStopWatch CGTimer;
  
  // SchurDiagTwoOperator<MobiusFermionD, LatticeFermionD> HermOpEO_outer(D_outer);

  // CGTimer.Start();
  // CG_outer(HermOpEO_outer, src_o_outer, result_o_outer);
  // CGTimer.Stop();

  // std::cout << GridLogMessage << "Total outer CG time : " << CGTimer.Elapsed()
  //           << std::endl;

  // CGTimer.Reset();

  //Solve for y using MADWF with internal preconditioning

  //typedef PauliVillarsSolverRBprec<LatticeFermionD, typename RunParamsOuter::SchurSolverType> PVtype;
  //PVtype PV_outer(SchurSolver_outer);

  typedef PauliVillarsSolverFourierAccel<LatticeFermionD, LatticeGaugeFieldD> PVtype;
  PVtype PV_outer(Umu, CG_outer);

  ConjugateGradient<LatticeFermionD> CG_inner(resid_inner, 10000, 0);

  CGincreaseTol update(CG_inner, resid_outer);

  SchurRedBlackDiagTwoSolve<LatticeFermionD> SchurSolver_inner(CG_inner);

  ZeroGuesser<LatticeFermion> Guess;
  MADWF<MobiusFermionD, ZMobiusFermionD, PVtype, SchurRedBlackDiagTwoSolve<LatticeFermionD>, ZeroGuesser<LatticeFermion> > madwf(D_outer, D_inner, PV_outer, SchurSolver_inner, Guess, resid_outer, 100, &update);
  
  // MADWF<MobiusFermionD, ZMobiusFermionD, PVtype, SchurRedBlackDiagTwoSolve<LatticeFermionD>, LinearFunction<LatticeFermion> > madwf(D_outer, D_inner, PV_outer, SchurSolver_inner, *guesserPt, resid_outer, 100, &update);


  LatticeFermionD result_MADWF(FGrid_outer);
  result_MADWF = Zero();

  CGTimer.Start();
  madwf(src4, result_MADWF);
  CGTimer.Stop();

  LatticeFermionD result_o_MADWF(FrbGrid_outer);
  pickCheckerboard(Odd, result_o_MADWF, result_MADWF);

  std::cout << GridLogMessage << "Total MADWF time : " << CGTimer.Elapsed()
            << std::endl;

  // LatticeFermionD diff = result_o_MADWF - result_o_outer;
  // std::cout <<GridLogMessage<< "Odd-parity MADWF result norm " << norm2(result_o_MADWF) 
  //       << " Regular result norm " << norm2(result_o_outer) 
  //       << " Norm of diff " << norm2(diff)<<std::endl;


  std::cout << GridLogMessage << "######## Dhop calls summary" << std::endl;
  D_outer.Report();

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_MADWFPrecCG_hpp_
