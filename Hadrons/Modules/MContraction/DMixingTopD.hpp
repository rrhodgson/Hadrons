/*
 * WeakEye3pt.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@cern.ch>
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
#ifndef Hadrons_MContraction_DMixingTopD_hpp_
#define Hadrons_MContraction_DMixingTopD_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DMixingTopD                                        *
 *                    (Fig. 4 (D) in arxiv:2504.16189)  
 *             qCL   ┌───┐                 qUR
 *         /----<----| r |-----------------<--------------\
 *        /          └───┘                                 \
 *       /          /     \          qLoop2 /->-\           \
 *   g5 *           \     /                /     \           * g5
 *       \           \->-/ qLoop1          \     /          /   
 *        \                                 ┌───┐          /       
 *         \-------------->-----------------| r'|---->----/  
 *                       qUL                └───┘  qCR
 *
 * four configurations for the two weak Hamiltonians M_r M_{r'}
 * (cf. Fig. 3 in arxiv:2504.16189)
 * rr'=11: tr()
 * rr'=12: tr()
 * rr'=21: tr()
 * rr'=22: tr()
 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class DMixingTopDPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMixingTopDPar,
                                    std::string,    qULeft,
                                    std::string,    qCLeft,
                                    std::string,    qURight,
                                    std::string,    qCRight,
                                    std::string,    qLoop1,
                                    std::string,    qLoop2,
                                    std::string,    output);
};

template <typename FImpl>
class TDMixingTopD: public Module<DMixingTopDPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        std::string,    parity,
                                        std::string,    rr);
    };
    typedef Correlator<Metadata> Result;
public:
    // constructor
    TDMixingTopD(const std::string name);
    // destructor
    virtual ~TDMixingTopD(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(DMixingTopD, TDMixingTopD<FIMPL>, MContraction);

/******************************************************************************
 *                        TDMixingTopD implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDMixingTopD<FImpl>::TDMixingTopD(const std::string name)
: Module<DMixingTopDPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDMixingTopD<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qULeft, 
	                           par().qCLeft,
	                           par().qURight,
	                           par().qCRight,
	                           par().qLoop1,
	                           par().qLoop2};

    return in;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopD<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopD<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;

    if (!par().output.empty())
        output.push_back(resultFilename(par().output));

    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopD<FImpl>::setup(void)
{
    envTmpLat(ComplexField, "corr");
    envTmpLat(PropagatorField, "parPlusR1L1");
    envTmpLat(PropagatorField, "parPlusR2L1");
    envTmpLat(PropagatorField, "parMinusR1L1");
    envTmpLat(PropagatorField, "parMinusR2L1");
    envTmpLat(PropagatorField, "parPlusR1L2");
    envTmpLat(PropagatorField, "parPlusR2L2");
    envTmpLat(PropagatorField, "parMinusR1L2");
    envTmpLat(PropagatorField, "parMinusR2L2");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopD<FImpl>::execute(void)
{
    LOG(Message) << "Computing D-meson mixing diagram, topology D" << std::endl;
    LOG(Message) << "qULeft  : " << par().qULeft << std::endl;
    LOG(Message) << "qCLeft  : " << par().qCLeft << std::endl;
    LOG(Message) << "qURight : " << par().qURight << std::endl;
    LOG(Message) << "qCRight : " << par().qCRight << std::endl;
    LOG(Message) << "qLoop1  : " << par().qLoop1 << std::endl;
    LOG(Message) << "qLoop2  : " << par().qLoop2 << std::endl;

    std::vector<Result> result;
    Result              r;

    auto                &qul  = envGet(PropagatorField, par().qULeft);
    auto                &qcl  = envGet(PropagatorField, par().qCLeft);
    auto                &qur  = envGet(PropagatorField, par().qURight);
    auto                &qcr  = envGet(PropagatorField, par().qCRight);
    auto                &ql1  = envGet(PropagatorField, par().qLoop1);
    auto                &ql2  = envGet(PropagatorField, par().qLoop2);

    Gamma               g5(Gamma::Algebra::Gamma5);
    Gamma               gVX(Gamma::Algebra::GammaX);
    Gamma               gVY(Gamma::Algebra::GammaY);
    Gamma               gVZ(Gamma::Algebra::GammaZ);
    Gamma               gVT(Gamma::Algebra::GammaT);
    Gamma               gAX(Gamma::Algebra::GammaXGamma5);
    Gamma               gAY(Gamma::Algebra::GammaYGamma5);
    Gamma               gAZ(Gamma::Algebra::GammaZGamma5);
    Gamma               gAT(Gamma::Algebra::GammaTGamma5);
    std::vector<Gamma>      GV = {gVX, gVY, gVZ, gVT};
    std::vector<Gamma>      GA = {gAX, gAY, gAZ, gAT};

    envGetTmp(ComplexField, corr);
    envGetTmp(PropagatorField, parPlusR1L1);
    envGetTmp(PropagatorField, parPlusR2L1);
    envGetTmp(PropagatorField, parMinusR1L1);
    envGetTmp(PropagatorField, parMinusR2L1);
    envGetTmp(PropagatorField, parPlusR1L2);
    envGetTmp(PropagatorField, parPlusR2L2);
    envGetTmp(PropagatorField, parMinusR1L2);
    envGetTmp(PropagatorField, parMinusR2L2);

    SlicedComplex buf;

    // rr'=11  (one single big trace)
    // for the moment, all gamma5 written explicitly - can simplyfy this
    // corr = tr(g5*g5*adj(qcl)*g5*G11*ql1*G12*qur*g5*g5*adj(qcr)*g5*G21*ql2*G22*qul); 
   
    // VV
    for (const auto &G: GV)
    {
        auto obj1 = G*ql1*G;
	parPlusR1L1 += obj1;
	obj1 = G*ql2*G;
	parPlusR1L2 += obj1;
        //auto obj2 = G*trace(ql1*G);
	//parPlusR2 += obj2;
    }
    // AA
    for (const auto &G: GA)
    {
        auto obj1 = G*ql1*G;
	parPlusR1L1 += obj1;
	obj1 = G*ql2*G;
	parPlusR1L2 += obj1;
    }
    parMinusR1L1 = parPlusR1L1 * g5;
    parMinusR1L2 = parPlusR1L2 * g5;

    // parity = + , rr = 11
    corr = trace(g5*g5*adj(qcl)*parPlusR1L1*qur*g5*g5*adj(qcr)*g5*parPlusR1L2*qul);
    sliceSum(corr, buf, Tp);
    r.corr.clear();
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        r.corr.push_back(TensorRemove(buf[t]));
    }
    r.info.parity = "+";
    r.info.rr     = "11";
    result.push_back(r);

    // parity = - , rr = 11
    corr = trace(g5*g5*adj(qcl)*parMinusR1L1*qur*g5*g5*adj(qcr)*g5*parMinusR1L2*qul);
    sliceSum(corr, buf, Tp);
    r.corr.clear();
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        r.corr.push_back(TensorRemove(buf[t]));
    }
    r.info.parity = "-";
    r.info.rr     = "11";
    result.push_back(r);

    // save result, and hand it to environment
    saveResult(par().output, "DMixingTopD", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DMixingTopD_hpp_
