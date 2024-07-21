/*
 * Z3.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
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

#ifndef Hadrons_MSource_Z3_hpp_
#define Hadrons_MSource_Z3_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Z_2 stochastic source
 -----------------------------
 * src_x = eta_x * theta(x_3 - tA) * theta(tB - x_3)
 
 the eta_x are independent uniform random numbers in {1,e^{i 2pi/3},e^{i 4pi/3}}
 
 * options:
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 
 */
 
/******************************************************************************
 *                          Z3 stochastic source                              *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class Z3Par: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Z3Par,
                                    unsigned int, tA,
                                    unsigned int, tB);
};

template <typename FImpl>
class TZ3: public Module<Z3Par>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TZ3(const std::string name);
    // destructor
    virtual ~TZ3(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool        hasT_{false};
    std::string tName_;
};

MODULE_REGISTER_TMP(Z3,       TZ3<FIMPL>,        MSource);
MODULE_REGISTER_TMP(ScalarZ3, TZ3<ScalarImplCR>, MSource);

/******************************************************************************
 *                       TZ3 template implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TZ3<FImpl>::TZ3(const std::string name)
: Module<Z3Par>(name)
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TZ3<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TZ3<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TZ3<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envTmpLat(LatticeComplex, "eta");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TZ3<FImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating Z_3 wall source at t= " << par().tA
                     << std::endl;
    }
    else
    {
        LOG(Message) << "Generating Z_3 band for " << par().tA << " <= t <= "
                     << par().tB << std::endl;
    }
    
    auto    &src = envGet(PropagatorField, getName());
    auto    &t   = envGet(Lattice<iScalar<vInteger>>, tName_);
    Complex i(0.0,1.0);

    if (!hasT_)
    {
        LatticeCoordinate(t, Tp);
        hasT_ = true;
    }
    envGetTmp(LatticeComplex, eta);

    unsigned int vol  = rng4d().Grid()->iSites()*rng4d().Grid()->oSites();
    std::vector<std::discrete_distribution<int32_t>> dist(vol,std::discrete_distribution<int32_t>{1,1,1});
    rng4d().fill(eta,dist);
    eta = real(eta);
    eta = exp((Real)(2*M_PI/3.)*i*eta);
    eta = where((t >= par().tA) and (t <= par().tB), eta, 0.*eta);
    src = 1.;
    src = src*eta;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Z3_hpp_
