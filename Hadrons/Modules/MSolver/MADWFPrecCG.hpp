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
                                    std::string , eigenPack);
};

template <typename FImplInner, typename FImplOuter, int nBasis>
class TMADWFPrecCG: public Module<MADWFPrecCGPar>
{
public:
    FERM_TYPE_ALIASES(FImplInner, Inner);
    FERM_TYPE_ALIASES(FImplOuter, Outer);
    SOLVER_TYPE_ALIASES(FImplOuter,);
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

MODULE_REGISTER_TMP(MADWFPrecCG, 
    ARG(TMADWFPrecCG<ZFIMPLF, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
// MODULE_REGISTER_TMP(ZMADWFPrecCG, 
//     ARG(TMADWFPrecCG<ZFIMPLD, ZFIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);

/******************************************************************************
 *                 TMADWFPrecCG implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
TMADWFPrecCG<FImplInner, FImplOuter, nBasis>
::TMADWFPrecCG(const std::string name)
: Module<MADWFPrecCGPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
std::vector<std::string> TMADWFPrecCG<FImplInner, FImplOuter, nBasis>
::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImplInner, typename FImplOuter, int nBasis>
std::vector<std::string> TMADWFPrecCG<FImplInner, FImplOuter, nBasis>
::getReference(void)
{
    std::vector<std::string> ref = {par().innerAction, par().outerAction};
    
    if (!par().eigenPack.empty())
    {
        ref.push_back(par().eigenPack);
    }
    
    return ref;
}

template <typename FImplInner, typename FImplOuter, int nBasis>
std::vector<std::string> TMADWFPrecCG<FImplInner, FImplOuter, nBasis>
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
template <typename FImplInner, typename FImplOuter, int nBasis>
void TMADWFPrecCG<FImplInner, FImplOuter, nBasis>
::setup(void)
{
    LOG(Message) << "Setting up Schur red-black preconditioned MADWF "
                 << "CG for inner/outer action '" << par().innerAction 
                 << "'/'" << par().outerAction << "', residual "
                 << par().residual << ", and maximum inner/outer iteration " 
                 << par().maxInnerIteration << "/" << par().maxOuterIteration
                 << std::endl;

    auto Ls_inner        = env().getObjectLs(par().innerAction);
    auto Ls_outer        = env().getObjectLs(par().outerAction);

    auto &imat     = envGet(FMatInner, par().innerAction);
    auto &omat     = envGet(FMatOuter, par().outerAction);

    MobiusFermionD D_outer  = *env().template getDerivedObject<FermionOperator<FImplOuter>,MobiusFermionD>(par().outerAction);
    ZMobiusFermionD D_inner = *env().template getDerivedObject<FermionOperator<FImplInner>,ZMobiusFermionD>(par().innerAction);

    auto guesserPt = makeGuesser<FImplInner, nBasis>(par().eigenPack);


    auto makeSolver = [D_outer, D_inner, &omat, guesserPt, Ls_inner, Ls_outer, this](bool subGuess) 
    {
        return [D_outer, D_inner, &omat, guesserPt, subGuess, Ls_inner, Ls_outer, this]
        (FermionFieldOuter &sol, const FermionFieldOuter &source) 
        {

            ConjugateGradient<LatticeFermionD> CG_outer(par().residual, 10000);
            ConjugateGradient<LatticeFermionD> CG_inner(par().residual, 10000, 0);

            // MobiusFermionD D_outer(omat.GaugeGrid(), omat.FermionGrid(), omat.FermionRedBlackGrid(), omat.GaugeGrid(), omat.GaugeRedBlackGrid(), omat.Mass, omat.par().M5, omat.par().par_b, omat.par().par_c);
            // ZMobiusFermionD D_inner(imat.U, imat.g5, imat.grb5, imat.g4, imat.grb4, imat.par().mass, imat.par().M5, imat.par().par_b, imat.par().par_c); 
                                //TODO:  ^ these member variables don't exist (find out how to accesss this info)

            LatticeFermionD src_outer(omat.FermionGrid());
            D_outer.ImportPhysicalFermionSource(source,src_outer); //applies D_- 
            //                             // TODO: ^ should src4 be source

            typedef PauliVillarsSolverFourierAccel<LatticeFermionD, LatticeGaugeFieldD> PVtype;
            PVtype PV_outer(omat.U, CG_outer);


            typedef SchurRedBlackDiagTwoSolve<LatticeFermionD> SchurSolverType;

            CGincreaseTol update(CG_inner, par().residual);

            SchurSolverType SchurSolver_inner(CG_inner);
            
            MADWF<MobiusFermionD, ZMobiusFermionD, PVtype, SchurSolverType, DeflatedGuesser<LatticeFermion> > 
                    madwf(D_outer, D_inner, PV_outer, SchurSolver_inner, *guesserPt, par().residual, 100, &update);
                    


            // typedef typename FermionFieldInner::vector_type VTypeInner;

            // SchurFMatInner simat(imat);
            // SchurFMatOuter somat(omat);
            // MixedPrecisionConjugateGradient<FermionFieldOuter, FermionFieldInner> 
            //     mpcg(par().residual, par().maxInnerIteration, 
            //          par().maxOuterIteration, 
            //          env().template getRbGrid<VTypeInner>(Ls),
            //          simat, somat);
            // OperatorFunctionWrapper<FermionFieldOuter> wmpcg(mpcg);
            // HADRONS_DEFAULT_SCHUR_SOLVE<FermionFieldOuter> schurSolver(wmpcg);
            // schurSolver.subtractGuess(subGuess);
            // schurSolver(omat, source, sol, *guesserPt);
        };
    };
    auto solver = makeSolver(false);
    envCreate(Solver, getName(), Ls_outer, solver, omat);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls_outer, solver_subtract, omat);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, int nBasis>
void TMADWFPrecCG<FImplInner, FImplOuter, nBasis>
::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_MADWFPrecCG_hpp_
