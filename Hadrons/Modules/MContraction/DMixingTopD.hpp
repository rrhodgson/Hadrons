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
                                        std::string,    rr,
                                        std::string,    parity,
					std::string,    eta_max);
    };
    typedef Correlator<Metadata, std::vector<ComplexD>> Result;
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
    // bespoke subcontractions
    virtual std::vector<SpinColourMatrixD> contract_D_half(const LatticeSpinColourMatrixD& prop_c, const LatticeSpinColourMatrixD& prop_u, const LatticeSpinColourMatrixD& loop);
    virtual std::vector<std::vector<ComplexD>> contract_D(const std::vector<SpinColourMatrixD>& half_if, const std::vector<SpinColourMatrixD>& half_fi);
    virtual std::pair<LatticeSpinColourMatrixD,LatticeSpinColourMatrixD> GH_VVAA_cap(const LatticeSpinColourMatrixD& prop);
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

template <typename FImpl>
std::vector<SpinColourMatrixD> TDMixingTopD<FImpl>::contract_D_half(const LatticeSpinColourMatrixD& prop_c, const LatticeSpinColourMatrixD& prop_u, const LatticeSpinColourMatrixD& loop) {
	Gamma G5(Gamma::Algebra::Gamma5);

	LatticePropagator tmp = G5*adj(prop_u)*G5 * loop * prop_c;
	std::vector<LatticePropagator::scalar_object> ret;
	sliceSum(tmp, ret, Tp);
	return ret;
};

template <typename FImpl>
std::vector<std::vector<ComplexD>> TDMixingTopD<FImpl>::contract_D(const std::vector<SpinColourMatrixD>& half_if, const std::vector<SpinColourMatrixD>& half_fi) {
	Gamma G5(Gamma::Algebra::Gamma5);
	Gamma GT(Gamma::Algebra::GammaT);

	// Kept general in case anyone ever wants to play with this
	Gamma Gsrc = G5;
	Gamma Gsnk = Gsrc; // no conj on final interpolator for D-Dbar mixing

	int Nt = half_if.size();

	std::vector<std::vector<ComplexD>> corr(Nt,std::vector<ComplexD>(Nt,0.));
	for (int t1=0; t1<Nt; t1++)
	{
	    for (int t2=0; t2<Nt; t2++)
	    {
                corr[t1][t2] = TensorRemove(trace( half_if[t1] * Gsrc * half_fi[t2] * Gsnk ));
	    }
	}
	
	return corr;
};

template <typename FImpl>
std::pair<LatticeSpinColourMatrixD,LatticeSpinColourMatrixD> TDMixingTopD<FImpl>::GH_VVAA_cap(const LatticeSpinColourMatrixD& prop) {
	GridBase* grid = prop.Grid();

	std::array<Gamma,8> GHs{Gamma(Gamma::Algebra::GammaX),
	              	Gamma(Gamma::Algebra::GammaY),
	              	Gamma(Gamma::Algebra::GammaZ),
	              	Gamma(Gamma::Algebra::GammaT),
	              	Gamma(Gamma::Algebra::GammaXGamma5),
	              	Gamma(Gamma::Algebra::GammaYGamma5),
	              	Gamma(Gamma::Algebra::GammaZGamma5),
	              	Gamma(Gamma::Algebra::GammaTGamma5)};

	SpinColourMatrixD spId = Zero();
	for (int s=0; s<4; s++)
	{
	    for (int c=0; c<3; c++)
	    {
		spId()(s,s)(c,c) = 1.;
	    }
	}

        LatticeSpinColourMatrixD GTrPropG_VVAA(grid); GTrPropG_VVAA = Zero();
	LatticeSpinColourMatrixD GPropG_VVAA(grid)  ; GPropG_VVAA   = Zero();
	for (int g=0; g<GHs.size(); g++) 
	{
		Gamma GH = GHs[g];
		GTrPropG_VVAA += spId * GH * trace(prop * GH);
		GPropG_VVAA   +=        GH *       prop * GH ;
	}
	return std::make_pair(GTrPropG_VVAA, GPropG_VVAA);
};

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

    if (par().qLoop1 != par().qLoop2)
    {
        HADRONS_ERROR(Argument, "Current implementation for identical loops only");
    }


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

    std::vector<Result>                 result;
    Result                              res;

    const int Nt{env().getDim(Tdir)};
    GridCartesian * grid = envGetGrid(FermionField);

    auto                &qul  = envGet(PropagatorField, par().qULeft);
    auto                &qcl  = envGet(PropagatorField, par().qCLeft);
    auto                &qur  = envGet(PropagatorField, par().qURight);
    auto                &qcr  = envGet(PropagatorField, par().qCRight);
    auto                &ql1  = envGet(std::vector<PropagatorField*>, par().qLoop1);
    auto                &ql2  = envGet(std::vector<PropagatorField*>, par().qLoop2);

    int Neta = ql1.size();

    std::vector<std::vector<std::vector<SpinColourMatrixD>>> half_lr(Neta, std::vector<std::vector<SpinColourMatrixD>>(2, std::vector<SpinColourMatrixD>(Nt)));
    std::vector<std::vector<std::vector<SpinColourMatrixD>>> half_rl(Neta, std::vector<std::vector<SpinColourMatrixD>>(2, std::vector<SpinColourMatrixD>(Nt)));



    // parity +, parity -
    std::vector<Gamma> parityG = {Gamma(Gamma::Algebra::Identity),Gamma(Gamma::Algebra::Gamma5)};
    std::vector<LatticeSpinColourMatrixD> GdsG_pp(2, grid);
    for (int p = 0; p < 2; p++) 
    {
        for (int i=0; i<Neta; i++) 
        {
            // here one has to add ql2 if one wants them to be allowed to be different
            auto tmp = GH_VVAA_cap(*ql1[i]);  	
            GdsG_pp[0] = tmp.first;  // r1 
            GdsG_pp[1] = tmp.second; // r2 
            for (int r = 0; r < 2; r++) 
	    {
	        half_lr[i][r] = contract_D_half(qcl, qur, GdsG_pp[i + Neta * (r-1)] * parityG.at(p));
	        half_rl[i][r] = contract_D_half(qcr, qul, GdsG_pp[i + Neta * (r-1)] * parityG.at(p));
	    }
	}
    
        for (int r = 0; r < 2; r++) 
        {
            for (int s = 0; s < 2; s++) 
	    {
                res.info.rr = std::to_string(r+1) + std::to_string(s+1);
                res.info.parity = (p == 0) ? "+" : "-";
		std::vector<std::vector<std::vector<std::vector<ComplexD>>>> buf(Neta, std::vector<std::vector<std::vector<ComplexD>>>(Neta, std::vector<std::vector<ComplexD>>(Nt, std::vector<ComplexD>(Nt))));
		for (int i=0; i<Neta; i++) 
		{
		    for (int j=0; j<Neta; j++) 
		    {
			buf[i][j] = contract_D(half_lr[i][r], half_rl[j][s]);
		    }
		}

		// Average noises up to imax (+ remove diagonal terms)
		for (int imax=1; imax<=Neta; imax++) 
		{
		    std::vector<std::vector<ComplexD>> tmp = std::vector<std::vector<ComplexD>>(Nt,std::vector<ComplexD>(Nt,0.));
		    for (int i=0; i<imax; i++) {
		        for (int j=0; j<imax; j++) {
			    if (i != j) {
			        const auto& c = buf[i][j];
				for (int t1=0; t1<Nt; t1++) 
				{
				    for (int t2=0; t2<Nt; t2++) 
				    {
				        tmp[t1][t2] += c[t1][t2];
				    }
				}
			    }
			}
		    }
		    if (imax > 1) 
		    {
		        for (int t1=0; t1<Nt; t1++) 
			{
			    for (int t2=0; t2<Nt; t2++) 
			    {
			        tmp[t1][t2] /= imax*(imax-1); 
			    }
			}
		    }
		    res.info.eta_max = std::to_string(imax);
                    res.corr.clear();
		    res.corr = tmp;
                    result.push_back(res);
		}
	    }
	}
    }

    // save result, and hand it to environment
    saveResult(par().output, "DMixingTopD", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;

}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DMixingTopD_hpp_
