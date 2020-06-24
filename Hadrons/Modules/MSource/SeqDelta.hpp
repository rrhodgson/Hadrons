/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSource/SeqDelta.hpp

Copyright (C) 2015-2019

Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef Hadrons_MSource_SeqDelta_hpp_
#define Hadrons_MSource_SeqDelta_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential delta
 -----------------------------
 * src_x = q_x * delta(x,position)
 
 * options:
 - q: input propagator (string)
 - position: 4-position in delta, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         SeqDelta                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SeqDeltaPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SeqDeltaPar,
                                    std::string,    q,
                                    std::string,    position);
};

template <typename FImpl>
class TSeqDelta: public Module<SeqDeltaPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSeqDelta(const std::string name);
    // destructor
    virtual ~TSeqDelta(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SeqDelta, TSeqDelta<FIMPL>, MSource);

/******************************************************************************
 *                         TSeqDelta implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSeqDelta<FImpl>::TSeqDelta(const std::string name)
: Module<SeqDeltaPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSeqDelta<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSeqDelta<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqDelta<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqDelta<FImpl>::execute(void)
{
    LOG(Message) << "Isolating position [" << par().position << "] in source " << par().q << std::endl;

    auto  &src = envGet(PropagatorField, getName());
    auto  &q   = envGet(PropagatorField, par().q);
    SitePropagator q_site;
    Coordinate site(strToVec<int>(par().position));
    src = Zero();

    peekSite(q_site, q , site);

    pokeSite(q_site, src , site);

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_SeqDelta_hpp_


// void peekSite(sobj &s,const Lattice<vobj> &l,const Coordinate &site);
// void pokeSite(const sobj &s,Lattice<vobj> &l,const Coordinate &site);