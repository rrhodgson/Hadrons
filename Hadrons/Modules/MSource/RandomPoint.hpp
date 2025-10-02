/*
 * RandomPoint.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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

#ifndef Hadrons_MSource_RandomPoint_hpp_
#define Hadrons_MSource_RandomPoint_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 RandomPoint source
 ------------
 * src_x = delta_x,position
 
 * options:
 - position: space-separated integer sequence (e.g. "0 1 1 0")
 
 */

/******************************************************************************
 *                                  TRandomPoint                                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class RandomPointPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RandomPointPar,
                                    std::string, position);
};

template <typename FImpl>
class TRandomPoint: public Module<RandomPointPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TRandomPoint(const std::string name);
    // destructor
    virtual ~TRandomPoint(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RandomPoint,       TRandomPoint<FIMPL>,        MSource);
MODULE_REGISTER_TMP(ScalarRandomPoint, TRandomPoint<ScalarImplCR>, MSource);

/******************************************************************************
 *                       TRandomPoint template implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TRandomPoint<FImpl>::TRandomPoint(const std::string name)
: Module<RandomPointPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TRandomPoint<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TRandomPoint<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(),getName()+"_pos",getName()+"_str"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandomPoint<FImpl>::setup(void)
{
    envCreateLat(PropagatorField , getName());
    envCreate(std::vector<int>, getName()+"_pos", 1, 0);
    envCreate(std::string     , getName()+"_str", 1, "");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandomPoint<FImpl>::execute(void)
{
    auto tmp = strToVec<std::string>(par().position);
    int Nd = tmp.size();

    if (Nd != env().getNd())
    {
        HADRONS_ERROR(Size, "position has " + std::to_string(Nd)
                      + " components (must have " + std::to_string(env().getNd()) + ")");
    }

    auto& rng = rngSerial();

    auto &position = envGet(std::vector<int>, getName()+"_pos");
    position.resize(Nd);
    for (int d=0; d<Nd; d++)
    {
        if (tmp[d] == "*") 
        {
            int NL = env().getDim(d);
            auto dist = std::uniform_int_distribution<uint32_t>(0,NL-1);
            fillScalar(position[d], dist, rng._generators[0]);
        } 
        else 
        {
            position[d] = stoi(tmp[d]);
        }
    }

    auto &str = envGet(std::string, getName()+"_str");
    str = "";
    for (int d=0; d<Nd; d++)
    {
        str += std::to_string(position[d]);
        if (d != Nd-1)
            str += " ";
    }

    LOG(Message) << "Creating random point source at position " << position << std::endl;

    auto             &src     = envGet(PropagatorField, getName());
    SitePropagator   id;
    
    id  = 1.;
    src = Zero();
    pokeSite(id, src, position);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_RandomPoint_hpp_
