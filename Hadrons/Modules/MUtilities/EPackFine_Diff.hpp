/*
 * EPackFine_Diff.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: fionnoh <fionnoh@pscomp.ie>
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
#ifndef Hadrons_MUtilities_EPackFine_Diff_hpp_
#define Hadrons_MUtilities_EPackFine_Diff_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Load eigen vectors/values package                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class EPackFine_DiffPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(EPackFine_DiffPar,
                                    std::string, epack1,
                                    std::string, epack2);
};

template <typename Pack>
class TEPackFine_Diff: public Module<EPackFine_DiffPar>
{
public:
    typedef typename Pack::Field   Field;
    typedef typename Pack::FieldIo FieldIo;
    typedef BaseEigenPack<Field>   BasePack;

public:
    // constructor
    TEPackFine_Diff(const std::string name);
    // destructor
    virtual ~TEPackFine_Diff(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(EPackFine_Diff, ARG(TEPackFine_Diff<FermionEigenPack<FIMPL>>), MUtilities);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(EPackFine_DiffF, ARG(TEPackFine_Diff<FermionEigenPack<FIMPLF>>), MUtilities);
MODULE_REGISTER_TMP(EPackFine_DiffIo32, ARG(TEPackFine_Diff<FermionEigenPack<FIMPL, FIMPLF>>), MUtilities);
MODULE_REGISTER_TMP(EPackFine_DiffFIo64, ARG(TEPackFine_Diff<FermionEigenPack<FIMPLF, FIMPL>>), MUtilities);
#endif

/******************************************************************************
 *                    TEPackFine_Diff implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack>
TEPackFine_Diff<Pack>::TEPackFine_Diff(const std::string name)
: Module<EPackFine_DiffPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack>
std::vector<std::string> TEPackFine_Diff<Pack>::getInput(void)
{
    std::vector<std::string> in = {par().epack1, par().epack2};
    
    return in;
}

template <typename Pack>
std::vector<std::string> TEPackFine_Diff<Pack>::getOutput(void)
{
    std::vector<std::string> out;
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack>
void TEPackFine_Diff<Pack>::setup(void)
{
    // GridBase  *grid, *gridIo = nullptr, *gridRb = nullptr;

    // grid   = getGrid<Field>(par().Ls);
    // gridRb = getGrid<Field>(par().redBlack, par().Ls);
    // if (typeHash<Field>() != typeHash<FieldIo>())
    // {
    //     gridIo = getGrid<FieldIo>(par().redBlack, par().Ls);
    // }
    // envCreateDerived(BasePack, Pack, getName(), par().Ls, par().size, gridRb, gridIo);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack>
void TEPackFine_Diff<Pack>::execute(void)
{
    auto &epack1 = envGet(BasePack, par().epack1);
    auto &epack2 = envGet(BasePack, par().epack2);

    if (epack1.evec.size() == epack2.evec.size()) {

        int size = epack1.evec.size();

        for (int i=0; i<size; i++) {
            LOG(Message) << "Eigenvector " << i << std::endl;

            Field tot = epack1.evec[i] + epack2.evec[i];
            Field diff = epack1.evec[i] - epack2.evec[i];

            Real tot_mod = norm2(tot);
            Real diff_mod = norm2(diff);
            Complex ip  = innerProduct(epack1.evec[i] , epack2.evec[i]);

            LOG(Message) << "||ev1||^2 = " << norm2(epack1.evec[i]) << std::endl;
            LOG(Message) << "||ev2||^2 = " << norm2(epack2.evec[i]) << std::endl;
            LOG(Message) << "||ev1+ev2||^2 = " << tot_mod << std::endl;
            LOG(Message) << "||ev1-ev2||^2 = " << diff_mod << std::endl;
            LOG(Message) << "ev1.ev2 = " << ip << "   mod = " << abs(ip) << "   arg = " << arg(ip) << std::endl;
            LOG(Message) << std::endl;
        }

    } else {
        LOG(Message) << "Size of epack1 (" << epack1.evec.size() << ") != epack2 (" << epack2.evec.size() << ")" << std::endl;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_EPackFine_Diff_hpp_
