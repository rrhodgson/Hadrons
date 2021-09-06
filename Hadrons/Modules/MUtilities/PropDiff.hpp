/*
 * PropDiff.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MUtilities_PropDiff_hpp_
#define Hadrons_MUtilities_PropDiff_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                PropDiff                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class PropDiffPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PropDiffPar,
                                    std::string, source1,
                                    std::string, source2);
};

template <typename FImpl>
class TPropDiff: public Module<PropDiffPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TPropDiff(const std::string name);
    // destructor
    virtual ~TPropDiff(void) {};
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

MODULE_REGISTER_TMP(PropDiff, TPropDiff<FIMPL>, MUtilities);

/******************************************************************************
 *                      TPropDiff implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPropDiff<FImpl>::TPropDiff(const std::string name)
: Module<PropDiffPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPropDiff<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source1, par().source2};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TPropDiff<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropDiff<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropDiff<FImpl>::execute(void)
{
    LOG(Message) << "Computing quark propagator difference '" << getName() << "' =  '" << par().source1 << "' - '" << par().source2 << "'"
                 << std::endl;
    

    auto &propdiff = envGet(PropagatorField, getName());
    auto &prop1    = envGet(PropagatorField, par().source1);
    auto &prop2    = envGet(PropagatorField, par().source2);

    propdiff = prop1 - prop2;

    LOG(Message) << "||" << par().source1 << "||^2 = " << norm2(prop1) << std::endl;
    LOG(Message) << "||" << par().source2 << "||^2 = "<< norm2(prop2) << std::endl;
    LOG(Message) << "||" << getName() << "||^2 = "<< norm2(propdiff) << std::endl;
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_PropDiff_hpp_
