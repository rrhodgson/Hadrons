/*
 * EPackDiff.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Guido Cossu <guido.cossu@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Nils Asmussen <n.asmussen@soton.ac.uk>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: pretidav <david.preti@csic.es>
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

#ifndef Hadrons_MUtilities_EPackDiff_hpp_
#define Hadrons_MUtilities_EPackDiff_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                EPackDiff                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class EPackDiffPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(EPackDiffPar,
                                    std::string, epack1,
                                    std::string, epack2);
};

template <typename Pack, typename FImpl, typename GImpl>
class TEPackDiff: public Module<EPackDiffPar>
{
public:
    typedef typename Pack::Field   Field;
    typedef typename Pack::FieldIo FieldIo;
    typedef BaseEigenPack<Field>   BasePack;
public:
    FERM_TYPE_ALIASES(FImpl,);
    // typedef typename Pack::PropagatorField               PropagatorField;
    // typedef typename Pack::Field               FermionField;
public:
    // constructor
    TEPackDiff(const std::string name);
    // destructor
    virtual ~TEPackDiff(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
};

// MODULE_REGISTER_TMP(EPackDiff, TEPackDiff<Pack, FImpl, GImpl>, MUtilities);

MODULE_REGISTER_TMP(EPackDiff, ARG(TEPackDiff<FermionEigenPack<FIMPL>, FIMPL, GIMPL>), MUtilities);
MODULE_REGISTER_TMP(EPackDiffF, ARG(TEPackDiff<FermionEigenPack<FIMPLF>, FIMPLF, GIMPLF>), MUtilities);

/******************************************************************************
 *                      TEPackDiff implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack, typename FImpl, typename GImpl>
TEPackDiff<Pack, FImpl, GImpl>::TEPackDiff(const std::string name)
: Module<EPackDiffPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack, typename FImpl, typename GImpl>
std::vector<std::string> TEPackDiff<Pack, FImpl, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().epack1, par().epack2};
    
    return in;
}

template <typename Pack, typename FImpl, typename GImpl>
std::vector<std::string> TEPackDiff<Pack, FImpl, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack, typename FImpl, typename GImpl>
void TEPackDiff<Pack, FImpl, GImpl>::setup(void)
{
    envCreateLat(Field, "EVecdiff");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack, typename FImpl, typename GImpl>
void TEPackDiff<Pack, FImpl, GImpl>::execute(void)
{
    LOG(Message) << "Computing epack difference '" << getName() << "' =  '" << par().epack1 << "' - '" << par().epack2 << "'"
                 << std::endl;

    auto &epack1 = envGetDerived(BasePack, Pack, par().epack1);
    auto &epack2 = envGetDerived(BasePack, Pack, par().epack2);

    auto &EVecDiff = envGet(Field, "EVecdiff");
    EVecDiff.reset(epack1.evec[0].Grid());

    if (epack1.eval.size() == epack2.eval.size()) {
        int N = epack1.eval.size();
        for (int i=0; i<N; i++) {
            LOG(Message) << "val1[" << i << "] = " << epack1.eval[i] << "   "
                         << "val2[" << i << "] = " << epack2.eval[i] << "   "
                         << "diff = "<< epack1.eval[i] - epack2.eval[i] << std::endl;

            EVecDiff = epack1.evec[i] - epack2.evec[i];
            LOG(Message) << "||vec1[" << i << "]||^2 = " << norm2(epack1.evec[i]) << "   "
                         << "||vec2[" << i << "]||^2 = " << norm2(epack2.evec[i]) <<"   "
                         << "||diff||^2 = "<< norm2(EVecDiff) << std::endl << std::endl;
        }

    } else {
        LOG(Message) << "EPack size doesn't match!" << std::endl;
        LOG(Message) << "epack1 = " << epack1.eval.size() << std::endl;
        LOG(Message) << "epack2 = " << epack2.eval.size() << std::endl;
    }
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_EPackDiff_hpp_
