/*
 * MADWFCG.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>
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
#ifndef Hadrons_MSolver_MADWFCG_hpp_
#define Hadrons_MSolver_MADWFCG_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MSolver/Guesser.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                              MADWF CG Solver                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class MADWFCGPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MADWFCGPar,
                                    std::string , innerAction,
                                    std::string , outerAction,
                                    std::string , cleanupAction,
                                    unsigned int, maxInnerIteration,
                                    unsigned int, maxOuterIteration,
                                    unsigned int, maxPVIteration,
                                    double      , innerResidual,
                                    double      , outerResidual,
                                    std::string , eigenPack);
};

template <typename FImplInner, typename FImplOuter, typename FImplCleanUp, int nBasis>
class TMADWFCG: public Module<MADWFCGPar>
{
public:
    FERM_TYPE_ALIASES(FImplInner, Inner);
    FERM_TYPE_ALIASES(FImplOuter, Outer);
    FERM_TYPE_ALIASES(FImplCleanUp, CleanUp);
    SOLVER_TYPE_ALIASES(FImplCleanUp,);
public:
    // constructor
    TMADWFCG(const std::string name);
    // destructor
    virtual ~TMADWFCG(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    struct CGincreaseTol;
};

MODULE_REGISTER_TMP(ZMADWFCG,              ARG(TMADWFCG<ZFIMPLD, FIMPLD, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
MODULE_REGISTER_TMP(ZMADWFCGMixedPrec,     ARG(TMADWFCG<ZFIMPLF, FIMPLD, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);

MODULE_REGISTER_TMP(ZMADWFCGMoreMixedPrec, ARG(TMADWFCG<ZFIMPLF, FIMPLF, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);

MODULE_REGISTER_TMP( MADWFCG,              ARG(TMADWFCG< FIMPLD, FIMPLD, FIMPLD, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);

/******************************************************************************
 *                        TMADWFCG implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, typename FImplCleanUp, int nBasis>
TMADWFCG<FImplInner, FImplOuter, FImplCleanUp, nBasis>
::TMADWFCG(const std::string name)
: Module<MADWFCGPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, typename FImplCleanUp, int nBasis>
std::vector<std::string> TMADWFCG<FImplInner, FImplOuter, FImplCleanUp, nBasis>
::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImplInner, typename FImplOuter, typename FImplCleanUp, int nBasis>
std::vector<std::string> TMADWFCG<FImplInner, FImplOuter, FImplCleanUp, nBasis>
::getReference(void)
{
    std::vector<std::string> ref = {par().innerAction, par().outerAction, par().cleanupAction};
    
    if (!par().eigenPack.empty())
    {
        ref.push_back(par().eigenPack);
    }
    
    return ref;
}

template <typename FImplInner, typename FImplOuter, typename FImplCleanUp, int nBasis>
std::vector<std::string> TMADWFCG<FImplInner, FImplOuter, FImplCleanUp, nBasis>
::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_subtract"};

    return out;
}


template <typename FImplInner, typename FImplOuter, typename FImplCleanUp, int nBasis>
struct TMADWFCG<FImplInner, FImplOuter, FImplCleanUp, nBasis>
::CGincreaseTol : public MADWFinnerIterCallbackBase {
    ConjugateGradient<FermionFieldInner> &cg_inner;  
    RealD outer_resid;

    CGincreaseTol(ConjugateGradient<FermionFieldInner> &cg_inner,
    RealD outer_resid): cg_inner(cg_inner), outer_resid(outer_resid){}

    void operator()(const RealD current_resid){
        LOG(Message) << "CGincreaseTol with current residual " << current_resid << " changing inner tolerance " << cg_inner.Tolerance << " -> ";
        while(cg_inner.Tolerance < current_resid) cg_inner.Tolerance *= 2;

        //cg_inner.Tolerance = outer_resid/current_resid;
        LOG(Message) << cg_inner.Tolerance << std::endl;
    }
};

// template <class Fieldi, class Fieldo,IfNotSame<Fieldi,Fieldo> X=0>
// inline void convert(const Fieldi &from,Fieldo &to) 
// {
//   precisionChange(to,from);
// }
// template <class Fieldi, class Fieldo,IfSame<Fieldi,Fieldo> X=0>
// inline void convert(const Fieldi &from,Fieldo &to) 
// {
//   to=from;
// }



// setup ///////////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, typename FImplCleanUp, int nBasis>
void TMADWFCG<FImplInner, FImplOuter, FImplCleanUp, nBasis>
::setup(void)
{
    LOG(Message) << "Setting up MADWF solver " << std::endl;
    LOG(Message) << "with inner/outer action  '"          << par().innerAction       << "'/'" << par().outerAction       << std::endl;
    LOG(Message) << "     inner/outer residual "          << par().innerResidual     <<  "/"  << par().outerResidual     << std::endl;
    LOG(Message) << "     maximum inner/outer/PV iterations " << par().maxInnerIteration <<  "/"  << par().maxOuterIteration <<  "/"  << par().maxPVIteration << std::endl;

    auto Ls_outer  = env().getObjectLs(par().outerAction);
    auto &omat     = envGet(FMatCleanUp, par().cleanupAction);
    auto guesserPt = makeGuesser<FImplInner, nBasis>(par().eigenPack);

    CayleyFermion5D<FImplOuter> &D_outer = envGetDerived(FMatOuter, CayleyFermion5D<FImplOuter>, par().outerAction);
    CayleyFermion5D<FImplInner> &D_inner = envGetDerived(FMatInner, CayleyFermion5D<FImplInner>, par().innerAction);
    CayleyFermion5D<FImplCleanUp> &D_cleanup = envGetDerived(FMatCleanUp, CayleyFermion5D<FImplCleanUp>, par().cleanupAction);

    auto makeSolver = [&D_outer, &D_cleanup, &D_inner, guesserPt, Ls_outer, this] (bool subGuess)
    {
        return [&D_outer, &D_cleanup, &D_inner, guesserPt, Ls_outer, subGuess, this]
        (FermionFieldCleanUp &sol, const FermionFieldCleanUp &source)
        {
            if (subGuess) {
                HADRONS_ERROR(Implementation, "MADWF solver with subtracted guess is not implemented!");
            }

            LOG(Message) << "0 " << std::endl;

            ConjugateGradient<FermionFieldOuter> CG_PV(par().outerResidual, par().maxPVIteration);
            HADRONS_DEFAULT_SCHUR_SOLVE<FermionFieldOuter> Schur_PV(CG_PV);
            typedef PauliVillarsSolverRBprec<FermionFieldOuter, HADRONS_DEFAULT_SCHUR_SOLVE<FermionFieldOuter>> PVtype;
            PVtype PV_outer(Schur_PV);

            LOG(Message) << "1 " << std::endl;

            ConjugateGradient<FermionFieldInner> CG_inner(par().innerResidual, par().maxInnerIteration, 0);
            HADRONS_DEFAULT_SCHUR_SOLVE<FermionFieldInner> SchurSolver_inner(CG_inner);

            CGincreaseTol update(CG_inner, par().outerResidual);

            LOG(Message) << "2 " << std::endl;

            FermionFieldOuter pre_sol(D_outer.FermionGrid());
            // FermionFieldOuter pre_sol(env().getGrid(Ls_outer));
            convert(sol, pre_sol);
            FermionFieldOuter pre_source(D_outer.FermionGrid());
            // FermionFieldOuter pre_source(env().getGrid(Ls_outer));
            convert(source, pre_source);

            LOG(Message) << "3 " << std::endl;

            MADWF<CayleyFermion5D<FImplOuter>, CayleyFermion5D<FImplInner>,
                  PVtype, HADRONS_DEFAULT_SCHUR_SOLVE<FermionFieldInner>, 
                  LinearFunction<FermionFieldInner> >  
                madwf(D_outer, D_inner,
                      PV_outer, SchurSolver_inner,
                      *guesserPt,
                      par().outerResidual, par().maxOuterIteration,
                      &update);

            LOG(Message) << "4 " << std::endl;

            madwf(pre_source, pre_sol);

            LOG(Message) << "5 " << std::endl;

            convert(pre_sol, sol);

            LOG(Message) << "6 " << std::endl;

            ConjugateGradient<FermionFieldCleanUp> CG_correction(par().outerResidual, 10000);
              HADRONS_DEFAULT_SCHUR_SOLVE<FermionFieldCleanUp> shur_correction(CG_correction);

            LOG(Message) << "7 " << std::endl;

              FermionFieldCleanUp src_e_outer(env().getRbGrid(Ls_outer));
              FermionFieldCleanUp src_o_outer(env().getRbGrid(Ls_outer));
              shur_correction.RedBlackSource(D_cleanup, source, src_e_outer, src_o_outer);

            LOG(Message) << "8 " << std::endl;
            
              FermionFieldCleanUp sol_o_outer(env().getRbGrid(Ls_outer));
              // sol_o_outer = Zero();
              pickCheckerboard(Odd, sol_o_outer, sol);
              FermionFieldCleanUp Moosol_o_outer(env().getRbGrid(Ls_outer));
              D_cleanup.Mooee(sol_o_outer,Moosol_o_outer);

            LOG(Message) << "9 " << std::endl;
            

            LOG(Message) << " ||src_o||^2 = " << norm2(src_o_outer) << std::endl;
            LOG(Message) << " ||sol_o||^2 = " << norm2(sol_o_outer) << std::endl;
            LOG(Message) << " ||Mooee sol_o||^2 = " << norm2(Moosol_o_outer) << std::endl;

              SchurDiagTwoOperator<CayleyFermion5D<FImplCleanUp>,FermionFieldCleanUp> HermOpEO_outer(D_cleanup);

              CG_correction(HermOpEO_outer, src_o_outer, Moosol_o_outer);

            LOG(Message) << " ||sol_o||^2 = " << norm2(sol_o_outer) << std::endl;

              shur_correction.RedBlackSolution(D_cleanup,Moosol_o_outer,src_e_outer,sol);
        };
    };
    auto solver = makeSolver(false);
    envCreate(Solver, getName(), Ls_outer, solver, omat);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls_outer, solver_subtract, omat);
}


// execution ///////////////////////////////////////////////////////////////////
template <typename FImplInner, typename FImplOuter, typename FImplCleanUp, int nBasis>
void TMADWFCG<FImplInner, FImplOuter, FImplCleanUp, nBasis>
::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_MADWFCG_hpp_
