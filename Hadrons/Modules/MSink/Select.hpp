
/*
 * Select.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk.com>
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

#ifndef Hadrons_MSink_Select_hpp_
#define Hadrons_MSink_Select_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                 Select                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSink)

class SelectPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SelectPar,
                                    std::vector<std::vector<int>>, points);
};

template <typename FImpl>
class TSelect: public Module<SelectPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SINK_TYPE_ALIASES();
public:
    // constructor
    TSelect(const std::string name);
    // destructor
    virtual ~TSelect(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Select,       TSelect<FIMPL>,        MSink);

/******************************************************************************
 *                          TSelect implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSelect<FImpl>::TSelect(const std::string name)
: Module<SelectPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSelect<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSelect<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSelect<FImpl>::setup(void)
{
    envCreate(SinkFn, getName(), 1, nullptr);
}


// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSelect<FImpl>::execute(void)
{
    LOG(Message) << "Setting up sink function to select points " << par().points << std::endl;
    
    auto sink = [this](const PropagatorField &field)
    {
        int Nt = env().getDim(3);

        SlicedPropagator res(Nt, Zero());
        SitePropagator   tmp;

        for (int i=0; i<par().points.size(); i++) {
            peekSite(tmp, field, par().points[i]);
            res[par().points[i][3]] += tmp;
        }
        
        return res;
    };
    envGet(SinkFn, getName()) = sink;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSink_Select_hpp_