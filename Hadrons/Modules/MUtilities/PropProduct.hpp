/*
 * PropProduct.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
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
#ifndef Hadrons_MUtilities_PropProduct_hpp_
#define Hadrons_MUtilities_PropProduct_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         PropProduct                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class PropProductPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PropProductPar,
                                    std::string, q,
                                    std::string, qBar);
};

template <typename FImpl>
class TPropProduct: public Module<PropProductPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TPropProduct(const std::string name);
    // destructor
    virtual ~TPropProduct(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(PropProduct, TPropProduct<FIMPL>, MUtilities);

/******************************************************************************
 *                 TPropProduct implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPropProduct<FImpl>::TPropProduct(const std::string name)
: Module<PropProductPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPropProduct<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q,par().qBar};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TPropProduct<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropProduct<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropProduct<FImpl>::execute(void)
{
    LOG(Message) << par().q    << std::endl;
    LOG(Message) << par().qBar << std::endl;

    auto &res  = envGet(PropagatorField, getName());
    auto &q    = envGet(PropagatorField, par().q);
    auto &qBar = envGet(PropagatorField, par().qBar);

    Gamma g5(Gamma::Algebra::Gamma5);

    res = q*g5*adj(qBar)*g5;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_PropProduct_hpp_
