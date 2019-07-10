//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBEBLocation.h"
#include <cmath>

registerMooseObject("PhaseFieldApp", GBEBLocation);

template <>
InputParameters
validParams<GBEBLocation>()
{
  InputParameters params = validParams<DerivativeParsedMaterialHelper>();
  params.addClassDescription("Calculate value of grain boundaries in a polycrystalline sample");
  params.addParam<int>("s",20,"Sharpness of Grain Boundary edge");
  params.addParam<int>("t",20,"Thickness of Grain Boundary");
  params.addRequiredCoupledVarWithAutoBuild(
      "v","var_name_base","op_num","Array of coupled variables");

  return params;
}

GBEBLocation::GBEBLocation(const InputParameters & parameters)
  : DerivativeParsedMaterialHelper(parameters), 
    _op_num(coupledComponents("v")),
    _vals(_op_num),
    _thickness(getParam<int>("t")),
    _sharpness(getParam<int>("s"))
{
  for (unsigned int i = 0; i < _op_num; ++i)
    {
      EBTerm _val(getVar("v",i)->name().c_str());
      _vals[i] = _val;
    }
  EBFunction _GB_location;
  _GB_location(_vals[0]) = 1.0 + pow(_vals[0],_thickness);
  for (unsigned int i = 1; i < _op_num; ++i)
    {
      std::vector<EBTerm> _input_variables(_vals.begin(),_vals.begin() + i);
      std::vector<EBTerm> _output_variables(_vals.begin(),_vals.begin() + i + 1);
      _GB_location(_output_variables) = _GB_location(_input_variables) + pow(_vals[i],_thickness);
    }
  _GB_location(_vals) = pow(2.0/(1.0 + _GB_location(_vals)),_sharpness);

  functionParse(_GB_location); 

}
