/*
 * Baryon.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <dc-erbe1@tesseract-login1.ib0.sgi.cluster.dirac.ed.ac.uk>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>
 * Author: ferben <ferben@debian.felix.com>
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

#ifndef Hadrons_MContraction_Baryon_Dump_hpp_
#define Hadrons_MContraction_Baryon_Dump_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/BaryonUtils.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               Baryon                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaAB;
typedef std::pair<GammaAB, GammaAB> GammaABPair;

class Baryon_DumpPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Baryon_DumpPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, q3,
                                    std::string, quarks,
                                    std::string, shuffle,
                                    int, parity,
                                    std::string, output);
};

template <typename FImpl>
class TBaryon_Dump: public Module<Baryon_DumpPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SINK_TYPE_ALIASES();
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, gammaA_left,
                                        Gamma::Algebra, gammaB_left,
                                        Gamma::Algebra, gammaA_right,
                                        Gamma::Algebra, gammaB_right,
                                        std::string, quarksR,
                                        std::string, quarksL,
                                        std::string, shuffle,
                                        int, parity);
    };
    typedef Correlator<Metadata> Result;
public:
    // constructor
    TBaryon_Dump(const std::string name);
    // destructor
    virtual ~TBaryon_Dump(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // Which gamma algebra was specified
    Gamma::Algebra  al;
};

MODULE_REGISTER_TMP(Baryon_Dump, ARG(TBaryon_Dump<FIMPL>), MContraction);

/******************************************************************************
 *                         TBaryon_Dump implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TBaryon_Dump<FImpl>::TBaryon_Dump(const std::string name)
: Module<Baryon_DumpPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TBaryon_Dump<FImpl>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().q3};
    
    return input;
}

template <typename FImpl>
std::vector<std::string> TBaryon_Dump<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}
template <typename FImpl>
std::vector<std::string> TBaryon_Dump<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;

    int L = env().getDim(0);
    for(int x=0; x<L; x++) {
    for(int y=0; y<L; y++) {
    for(int z=0; z<L; z++) {
        resultFilename(par().output+"_x"+std::to_string(x)+"_y"+std::to_string(y)+"_z"+std::to_string(z));
    }}}
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBaryon_Dump<FImpl>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBaryon_Dump<FImpl>::execute(void)
{
    // Check shuffle is a permutation of "123"
    assert(par().shuffle.size()==3 && "shuffle parameter must be 3 characters long");
    std::string shuffle_tmp = par().shuffle;
    std::sort(shuffle_tmp.begin(), shuffle_tmp.end());
    assert(shuffle_tmp == "123" && "shuffle parameter must be a permulation of 123");

    std::vector<int> shuffle = { std::stoi( par().shuffle.substr(0,1) ) -1,
                                 std::stoi( par().shuffle.substr(1,1) ) -1,
                                 std::stoi( par().shuffle.substr(2,1) ) -1 };


    assert(par().quarks.size()==3 && "quark-structure must consist of 3 quarks");
    std::string quarksR = par().quarks;
    std::string quarksL = quarksR;

    // Shuffle quark flavours
    quarksL[0] = quarksR[shuffle[0]];
    quarksL[1] = quarksR[shuffle[1]];
    quarksL[2] = quarksR[shuffle[2]];

    std::vector<std::string> props(3, "");
    props[0] = par().q1;
    props[1] = par().q2;
    props[2] = par().q3;

    // Shuffle propagators
    std::vector<std::string> propsL(props);
    propsL[0] = props[shuffle[0]];
    propsL[1] = props[shuffle[1]];
    propsL[2] = props[shuffle[2]];

    assert(par().parity == 1 || par().parity == -1 && "parity must be 1 or -1");



    LOG(Message) << "Computing baryon_Dump contractions '" << getName() << "'" << std::endl;
    LOG(Message) << "  using shuffle " << shuffle << " and parity " << par().parity << std::endl;
    LOG(Message) << "  using quarksL (" << quarksL << ") with left propagators (" << propsL[0] << ", " << propsL[1] << ", and " << propsL[2] << ")" << std::endl;
    LOG(Message) << "  using quarksR (" << quarksR << ") ";
    LOG(Message) << "with simultaneous sink " << std::endl;
    
    LOG(Message) << "    with (Gamma^A,Gamma^B)_left = ( " << Gamma::Algebra::Identity << " , " << Gamma::Algebra::MinusGammaZGamma5 << "') and (Gamma^A,Gamma^B)_right = ( " << Gamma::Algebra::Identity << " , " << Gamma::Algebra::MinusGammaZGamma5 << ")" << std::endl;
    
    envGetTmp(LatticeComplex, c);
    int nt = env().getDim(Tp);

    Result              r;
    r.info.parity  = par().parity;

    r.info.quarksR  = quarksR;
    r.info.quarksL  = quarksL;
    r.info.shuffle = par().shuffle;
        
    bool wick_contractions[6];
    BaryonUtils<FIMPL>::Wick_Contractions(quarksL,quarksR,wick_contractions);
    
    PropagatorField &q1  = envGet(PropagatorField, propsL[0]);
    PropagatorField &q2  = envGet(PropagatorField, propsL[1]);
    PropagatorField &q3  = envGet(PropagatorField, propsL[2]);

    std::vector<TComplex> buf;

    r.info.gammaA_left = Gamma::Algebra::Identity;
    r.info.gammaB_left = Gamma::Algebra::MinusGammaZGamma5;
    r.info.gammaA_right = Gamma::Algebra::Identity;
    r.info.gammaB_right = Gamma::Algebra::MinusGammaZGamma5;

    Gamma gAl(Gamma::Algebra::Identity);
    Gamma gBl(Gamma::Algebra::MinusGammaZGamma5);
    Gamma gAr(Gamma::Algebra::Identity);
    Gamma gBr(Gamma::Algebra::MinusGammaZGamma5);

    c=Zero();
    BaryonUtils<FIMPL>::ContractBaryons(q1,q2,q3,
                                        gAl,gBl,gAr,gBr,
                                        wick_contractions,
                                        par().parity,
                                        c);

    int L = env().getDim(0);
    int T = env().getDim(3);

    for(int x=0; x<L; x++) {
    for(int y=0; y<L; y++) {
    for(int z=0; z<L; z++) {
        std::vector<Result> result;

        std::cout << "(x,y,z) = (" << x << ", " << y << ", " << z << ")" << std::endl; 

        // SinkFnScalar &sink = envGet(SinkFnScalar, par().sinkq1);
        Coordinate coor(4);
        coor[0]=x;
        coor[1]=y;
        coor[2]=z;

        SinkFnScalar sink = [&coor,T](const PropagatorFieldScalar &field)
        {
            SlicedPropagatorScalar res(T);

            for (int t=0; t<T; t++) {
                coor[3]=t;
                peekSite(res[t],field,coor);
            }             
            return res;
        };

        buf = sink(c);

        r.corr.clear();
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            r.corr.push_back(TensorRemove(buf[t]));
        }
        result.push_back(r);

        saveResult(par().output+"_x"+std::to_string(x)+"_y"+std::to_string(y)+"_z"+std::to_string(z), "baryon", result);

    }}}

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Baryon_Dump_hpp_
