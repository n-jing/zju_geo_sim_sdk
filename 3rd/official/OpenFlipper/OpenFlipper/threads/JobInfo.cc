//=============================================================================
//
//                               OpenFlipper
//        Copyright (C) 2008 by Computer Graphics Group, RWTH Aachen
//                           www.openflipper.org
//
//-----------------------------------------------------------------------------
//
//                                License
//
//  OpenFlipper is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  OpenFlipper is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with OpenFlipper.  If not, see <http://www.gnu.org/licenses/>.
//
//-----------------------------------------------------------------------------
//
//   $Revision: 8597 $
//   $Author: kremer $
//   $Date: 2010-02-23 15:04:40 +0100 (Di, 23. Feb 2010) $
//
//=============================================================================


#include "JobInfo.hh"
#include <iostream>


JobInfo::JobInfo() :
  id(""),
  description(""),
  currentStep(0),
  minSteps(0),
  maxSteps(100),
  blocking(false),
  blockingWidget(0)
{
}
