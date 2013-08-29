// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Jens Allmer $
// $Authors: Jens Allmer $
// --------------------------------------------------------------------------

#ifndef OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_MATERCREATOR_H
#define OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_MATERCREATOR_H

#include <OpenMS/config.h>
#include "OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/Mater.h"
#include "OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/DefaultMater.h"
#include "OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/RandomMater.h"
#include "OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/SimpleMater.h"
#include "OpenMS/ANALYSIS/DENOVO/MSNOVOGEN/ZipMater.h"

namespace OpenMS
{
  class OPENMS_DLLAPI MaterCreator
  {
private:
	MaterCreator();
	MaterCreator(const MaterCreator& other);
	MaterCreator & operator=(const MaterCreator &other);

public:
	  static boost::shared_ptr<Mater> getInstance(const String & mater, double precursorMass, double precursorMassTolerance, std::vector<const Residue*> aaList)
	  {
		  if(mater == "RandomMater")
		    return boost::shared_ptr<RandomMater>(new RandomMater(precursorMass,precursorMassTolerance,aaList));
		  else {
			if(mater == "SimpleMater")
			  return boost::shared_ptr<SimpleMater>(new SimpleMater(precursorMass,precursorMassTolerance,aaList));
			else {
			  if(mater == "ZipMater")
				return boost::shared_ptr<ZipMater>(new ZipMater(precursorMass,precursorMassTolerance,aaList));
			  else
			    return boost::shared_ptr<DefaultMater>(new DefaultMater(precursorMass,precursorMassTolerance,aaList));
			}
		  }
	  }
  };
} // namespace

#endif // OPENMS__ANALYSIS_DENOVO_MSNOVOGEN_MATERCREATOR_H
