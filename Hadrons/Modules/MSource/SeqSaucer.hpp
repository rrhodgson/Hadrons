/*
 * SeqSaucer.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: fionnoh <fionnoh@gmail.com>
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

#ifndef Hadrons_MSource_SeqSaucer_hpp_
#define Hadrons_MSource_SeqSaucer_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential source
 -----------------------------
 * src_x = q_x * theta(x_3 - tA) * theta(tB - x_3) * gamma * exp(i x.mom)
 
 * options:
 - q: input propagator (string)
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - gamma: gamma product to insert (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         Sequential gamma source                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SeqSaucerPar: Serializable
{
public:
	GRID_SERIALIZABLE_CLASS_MEMBERS(SeqSaucerPar,
									std::string,    q,
									std::string,    qLoop,
									unsigned int,   tA,
									unsigned int,   tB,
									Gamma::Algebra, gamma);
};

template <typename FImpl>
class TSeqSaucer: public Module<SeqSaucerPar>
{
public:
	FERM_TYPE_ALIASES(FImpl,);
public:
	// constructor
	TSeqSaucer(const std::string name);
	// destructor
	virtual ~TSeqSaucer(void) {};
	// dependency relation
	virtual std::vector<std::string> getInput(void);
	virtual std::vector<std::string> getOutput(void);
protected:
	// setup
	virtual void setup(void);
	// execution
	virtual void execute(void);
private:
	void makeSource(PropagatorField &src, const PropagatorField &q);
private:
	bool        hasPhase_{false};
	std::string momphName_, tName_;
};

MODULE_REGISTER_TMP(SeqSaucer, TSeqSaucer<FIMPL>, MSource);
MODULE_REGISTER_TMP(ZSeqSaucer, TSeqSaucer<ZFIMPL>, MSource);

/******************************************************************************
 *                         TSeqSaucer implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSeqSaucer<FImpl>::TSeqSaucer(const std::string name)
: Module<SeqSaucerPar>(name)
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSeqSaucer<FImpl>::getInput(void)
{
	std::vector<std::string> in = {par().q, par().qLoop};
	
	return in;
}

template <typename FImpl>
std::vector<std::string> TSeqSaucer<FImpl>::getOutput(void)
{
	std::vector<std::string> out = {getName()};
	
	return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqSaucer<FImpl>::setup(void)
{
	envCreateLat(PropagatorField, getName());

    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqSaucer<FImpl>::execute(void)
{
	if (par().tA == par().tB)
	{
		LOG(Message) << "Generating " << par().gamma
					 << " sequential source(s) at t= " << par().tA << std::endl;
	} else {
		LOG(Message) << "Generating " << par().gamma
					 << " sequential source(s) for "
					 << par().tA << " <= t <= " << par().tB << std::endl;
	}

	auto  &src     = envGet(PropagatorField, getName()); 
	auto  &q       = envGet(PropagatorField, par().q);
	auto  &qLoop   = envGet(PropagatorField, par().qLoop);

	LOG(Message) << "Using initial propagator '" << par().q << "'" << std::endl;
	LOG(Message) << "and loop propagator '" << par().qLoop << "'" << std::endl;
	
  	GridBase *grid = q.Grid();

  	// std::cout << q << std::endl << std::endl << std::endl << std::endl;
  	// std::cout << "regular prop" << std::endl;

	auto  &t   = envGet(Lattice<iScalar<vInteger>>, tName_);
    LatticeCoordinate(t, Tp);

	q = where((t >= par().tA) and (t <= par().tB), q, 0.*q);

  	std::cout << q << std::endl << std::endl << std::endl << std::endl;
  	std::cout << "masked prop" << std::endl;

	Gamma g(par().gamma);

	autoView( vsrc , src, CpuWrite);
	autoView( vq , q, CpuRead);
	autoView( vqLoop , qLoop, CpuRead);
	thread_for(ss,grid->oSites(),{
		auto Dsrc = vsrc[ss];
		auto Dq = vq[ss];
		auto DqLoop = vqLoop[ss];

		// Dsrc = Zero();
		
		// Do computation

		// Dsrc = g * DqLoop * g * Dq;
		// Dsrc = g * Dq;

		// for (int a=0; a<3; a++){
		// for (int b=0; b<3; b++){
		// // for (int i=0; i<3; i++){

		// 	for (int alpha=0; alpha<4; alpha++) {
		// 	for (int beta=0; beta<4; beta++) {
		// 		Dsrc()(alpha,beta)(a,b) = Dq()(alpha,beta)(a,b);
		// 	}
		// 	}

		// }}
		Dsrc() = Dq();
	}  );//end loop over lattice sites


  	std::cout << "src = " << src << std::endl << std::endl << std::endl << std::endl;
  	std::cout << "seq prop" << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_SeqSaucer_hpp_
