/*
 * Epack_Compress.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MUtilities_Epack_Compress_hpp_
#define Hadrons_MUtilities_Epack_Compress_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Local coherence Lanczos eigensolver                     *
 *****************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class Epack_CompressPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Epack_CompressPar,
                                    std::string,   action,
                                    std::string,   blockSize,
                                    int,           fineSize,
                                    int,           coarseSize,
                                    std::string,   output,
                                    bool,          multiFile,
                                    std::string,   epack,
                                    int,           size);
};

template <typename FImpl, int nBasis, typename FImplIo = FImpl>
class TEpack_Compress: public Module<Epack_CompressPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef LocalCoherenceLanczos<typename FImpl::SiteSpinor, 
                                  typename FImpl::SiteComplex, 
                                  nBasis>                  LCL;
    typedef BaseFermionEigenPack<FImpl>                    BasePack;
    typedef CoarseFermionEigenPack<FImpl, nBasis, FImplIo> CoarsePack;
    typedef typename CoarsePack::Field                     Field;
    typedef typename CoarsePack::FieldIo                   FieldIo;
    typedef typename CoarsePack::CoarseField               CoarseField;
    typedef typename CoarsePack::CoarseFieldIo             CoarseFieldIo;
    typedef HADRONS_DEFAULT_SCHUR_OP<FMat, FermionField>   SchurFMat;

    typedef FermionEigenPack<FImpl>                        FinePack;
public:
    // constructor
    TEpack_Compress(const std::string name);
    // destructor
    virtual ~TEpack_Compress(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Epack_Compress, ARG(TEpack_Compress<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MUtilities);
MODULE_REGISTER_TMP(Epack_Compress30, ARG(TEpack_Compress<FIMPL, 30>), MUtilities);
MODULE_REGISTER_TMP(Epack_Compress150, ARG(TEpack_Compress<FIMPL, 150>), MUtilities);
MODULE_REGISTER_TMP(Epack_Compress400, ARG(TEpack_Compress<FIMPL, 400>), MUtilities);

MODULE_REGISTER_TMP(Epack_CompressF, ARG(TEpack_Compress<FIMPLF, HADRONS_DEFAULT_LANCZOS_NBASIS>), MUtilities);
MODULE_REGISTER_TMP(Epack_CompressF30, ARG(TEpack_Compress<FIMPLF, 30>), MUtilities);
MODULE_REGISTER_TMP(Epack_CompressF150, ARG(TEpack_Compress<FIMPLF, 150>), MUtilities);
MODULE_REGISTER_TMP(Epack_CompressF400, ARG(TEpack_Compress<FIMPLF, 400>), MUtilities);


/******************************************************************************
 *                 TEpack_Compress implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
TEpack_Compress<FImpl, nBasis, FImplIo>::TEpack_Compress(const std::string name)
: Module<Epack_CompressPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
std::vector<std::string> TEpack_Compress<FImpl, nBasis, FImplIo>::getInput(void)
{
    std::vector<std::string> in = {par().action, par().epack};
    
    return in;
}

template <typename FImpl, int nBasis, typename FImplIo>
std::vector<std::string> TEpack_Compress<FImpl, nBasis, FImplIo>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
void TEpack_Compress<FImpl, nBasis, FImplIo>::setup(void)
{
    LOG(Message) << "Setting up local coherence Lanczos eigensolver for"
                 << " action '" << par().action << "' (" << nBasis
                 << " eigenvectors)..." << std::endl;
    
    unsigned int Ls        = env().getObjectLs(par().action);
    auto         blockSize = strToVec<int>(par().blockSize);
    GridBase     *gridIo = nullptr, *gridCoarseIo = nullptr;

    auto cg  = envGetCoarseGrid(CoarseField, blockSize, Ls);

    LOG(Message) << "Coarse grid: " << cg->GlobalDimensions() << std::endl;
    envCreate(CoarsePack, getName(), Ls,
                     par().fineSize, par().coarseSize, envGetRbGrid(Field, Ls), cg,
                     gridIo, gridCoarseIo);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
void TEpack_Compress<FImpl, nBasis, FImplIo>::execute(void)
{
    auto &coarsePack   = envGet(CoarsePack, getName());
    auto &finePack     = envGet(BasePack, par().epack);

    int sizeFine = par().fineSize;
    int sizeCoarse = par().coarseSize;

    coarsePack.record = finePack.record;

    LOG(Message) << "Taking the first " << sizeFine << " fine vectors" << std::endl;
    for (int i=0; i<sizeFine; i++) {
        coarsePack.eval[i] = finePack.eval[i];
        coarsePack.evec[i] = finePack.evec[i];
    }

    if (!par().output.empty())
    {
        std::cout << "Write " << sizeFine << " fine vectors" << std::endl;
        coarsePack.writeFine(par().output, par().multiFile, vm().getTrajectory());
    }

    Lattice<typename FImpl::SiteComplex> InnerProd(( envGetCoarseGrid(CoarseField, strToVec<int>(par().blockSize), env().getObjectLs(par().action)) )); 
    LOG(Message) <<" Block Gramm-Schmidt pass 1"<<std::endl;
    blockOrthonormalize(InnerProd,coarsePack.evec);
    LOG(Message) <<" Block Gramm-Schmidt pass 2"<<std::endl;
    blockOrthonormalize(InnerProd,coarsePack.evec);

    LOG(Message) << "Projecting to coarse basis" << std::endl;
    for (int i=0; i<finePack.evec.size(); i++) {
        LOG(Message) << "evec " << i << std::endl;
        blockProject(coarsePack.evecCoarse[i], finePack.evec[i], coarsePack.evec);
        coarsePack.evalCoarse[i] = finePack.eval[i];
    }

    if (!par().output.empty())
    {
        std::cout << "Write " << sizeCoarse << " coarse vectors" << std::endl;
        coarsePack.writeCoarse(par().output, par().multiFile, vm().getTrajectory());
    }


}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_Epack_Compress_hpp_
