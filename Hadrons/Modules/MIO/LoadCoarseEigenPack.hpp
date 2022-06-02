/*
 * LoadCoarseEigenPack.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MIO_LoadCoarseEigenPack_hpp_
#define Hadrons_MIO_LoadCoarseEigenPack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Load local coherence eigen vectors/values package             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadCoarseEigenPackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadCoarseEigenPackPar,
                                    std::string, filestem,
                                    bool,         multiFile,
                                    unsigned int, sizeFine,
                                    unsigned int, sizeCoarse,
                                    unsigned int, Ls,
                                    std::vector<int>, blockSize);
};

template <typename Pack>
class TLoadCoarseEigenPack: public Module<LoadCoarseEigenPackPar>
{
public:
    typedef typename Pack::Field                Field;
    typedef typename Pack::FieldIo              FieldIo;
    typedef typename Pack::CoarseField          CoarseField;
    typedef typename Pack::CoarseFieldIo        CoarseFieldIo;
    typedef CoarseEigenPack<Field, CoarseField, FieldIo, CoarseFieldIo> BasePack;
    template <typename vtype> 
    using iImplScalar = iScalar<iScalar<iScalar<vtype>>>;
    typedef iImplScalar<typename Pack::Field::vector_type> SiteComplex;
public:
    // constructor
    TLoadCoarseEigenPack(const std::string name);
    // destructor
    virtual ~TLoadCoarseEigenPack(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack, 
                    ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>>), MIO);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(LoadCoarseFermionEigenPackF, 
                    ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPLF, HADRONS_DEFAULT_LANCZOS_NBASIS>>), MIO);
MODULE_REGISTER_TMP(LoadCoarseFermionEigenPackIo32, 
                    ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS, FIMPLF>>), MIO);
#endif

MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack150, 
                    ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL, 150>>), MIO);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack150F, 
                    ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPLF, 150>>), MIO);
MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack150Io32, 
                    ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL, 150, FIMPLF>>), MIO);
#endif

MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack400, 
                    ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL, 400>>), MIO);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack400F, 
                    ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPLF, 400>>), MIO);
MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack400Io32, 
                    ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL, 400, FIMPLF>>), MIO);
#endif

/******************************************************************************
 *                 TLoadCoarseEigenPack implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack>
TLoadCoarseEigenPack<Pack>::TLoadCoarseEigenPack(const std::string name)
: Module<LoadCoarseEigenPackPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack>
std::vector<std::string> TLoadCoarseEigenPack<Pack>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Pack>
std::vector<std::string> TLoadCoarseEigenPack<Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack>
void TLoadCoarseEigenPack<Pack>::setup(void)
{
    GridBase     *gridIo = nullptr, *gridCoarseIo = nullptr;

    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = envGetRbGrid(FieldIo, par().Ls);
    }
    if (typeHash<CoarseField>() != typeHash<CoarseFieldIo>())
    {
        gridCoarseIo = envGetCoarseGrid(CoarseFieldIo, par().blockSize, par().Ls);
    }
    envCreateDerived(BasePack, Pack, getName(), par().Ls, par().sizeFine,
                     par().sizeCoarse, envGetRbGrid(Field, par().Ls), 
                     envGetCoarseGrid(CoarseField, par().blockSize, par().Ls),
                     gridIo, gridCoarseIo);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack>
void TLoadCoarseEigenPack<Pack>::execute(void)
{
    auto                 cg     = envGetCoarseGrid(CoarseField, par().blockSize, par().Ls);
    auto                 &epack = envGetDerived(BasePack, Pack, getName());
    Lattice<SiteComplex> dummy(cg);

    epack.read(par().filestem, par().multiFile, vm().getTrajectory());
    LOG(Message) << "Block Gramm-Schmidt pass 1"<< std::endl;
    blockOrthogonalise(dummy, epack.evec);
    LOG(Message) << "Block Gramm-Schmidt pass 2"<< std::endl;
    blockOrthogonalise(dummy, epack.evec);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadCoarseEigenPack_hpp_
