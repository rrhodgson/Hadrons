/*
 * BoolTest.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
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
#ifndef Hadrons_MUtilities_BoolTest_hpp_
#define Hadrons_MUtilities_BoolTest_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/OptionalBool.cpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                          Optional Bool Test                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class BoolTestPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BoolTestPar,
                                    bool, b1,
                                    optionalBool, b2);
};

template <typename FieldIn, typename FieldOut>
class TBoolTest: public Module<BoolTestPar>
{
public:
    // constructor
    TBoolTest(const std::string name);
    // destructor
    virtual ~TBoolTest(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(BoolTest, 
                    ARG(TBoolTest<GIMPLD::GaugeField, GIMPLF::GaugeField>),
                    MUtilities);

/******************************************************************************
 *                     TBoolTest implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FieldIn, typename FieldOut>
TBoolTest<FieldIn, FieldOut>::TBoolTest(const std::string name)
: Module<BoolTestPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FieldIn, typename FieldOut>
std::vector<std::string> TBoolTest<FieldIn, FieldOut>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FieldIn, typename FieldOut>
std::vector<std::string> TBoolTest<FieldIn, FieldOut>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FieldIn, typename FieldOut>
void TBoolTest<FieldIn, FieldOut>::setup(void)
{
    
}

#define show(var) (var ? "true" : "false")

// execution ///////////////////////////////////////////////////////////////////
template <typename FieldIn, typename FieldOut>
void TBoolTest<FieldIn, FieldOut>::execute(void)
{
    try{
        LOG(Message) << "b1 (bool)         = " << show(par().b1) << std::endl;
    } catch ( InvalidOptionalException& e ) {
        std::cout << std::endl;
        LOG(Message) << "b1 - invalid access" << std::endl;
    }
    try{
        LOG(Message) << "b2 (optionalBool) = " << show(par().b2) << std::endl;
    } catch ( InvalidOptionalException& e ) {
        std::cout << std::endl;
        LOG(Message) << "b2 - invalid access" << std::endl;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_BoolTest_hpp_
