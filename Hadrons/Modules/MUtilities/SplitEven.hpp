/*
 * SplitEven.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MUtilities_SplitEven_hpp_
#define Hadrons_MUtilities_SplitEven_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SplitEven                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class SplitEvenPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SplitEvenPar,
                                    std::string, q1,
                                    std::string, q2);
};

template <typename FImpl>
class TSplitEven: public Module<SplitEvenPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSplitEven(const std::string name);
    // destructor
    virtual ~TSplitEven(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SplitEven, TSplitEven<FIMPL>, MUtilities);

/******************************************************************************
 *                 TSplitEven implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSplitEven<FImpl>::TSplitEven(const std::string name)
: Module<SplitEvenPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSplitEven<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q1,par().q2};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSplitEven<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSplitEven<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSplitEven<FImpl>::execute(void)
{
    LOG(Message) << par().q1 << std::endl;
    LOG(Message) << par().q2 << std::endl;

    auto &res  = envGet(PropagatorField, getName());
    auto &q1   = envGet(PropagatorField, par().q1);
    auto &q2   = envGet(PropagatorField, par().q2);

    Gamma g5(Gamma::Algebra::Gamma5);

    res = q1 * g5 * adj(q2) * g5;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_SplitEven_hpp_
