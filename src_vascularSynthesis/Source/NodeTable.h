/*=========================================================================
 
	Program: VascuSynth
	Module: $RCSfile: NodeTable.h,v $
	Language: C++
	Date: $Date: 2011/02/08 10:43:00 $
	Version: $Revision: 1.0 $

Copyright (c) 2011 Medical Imaging Analysis Lab, Simon Fraser University, 
British Columbia, Canada.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

* The name of the Insight Consortium, nor the names of any consortium members,
nor of any contributors, may be used to endorse or promote products derived
from this software without specific prior written permission.

* Modified source versions must be plainly marked as such, and must not be
misrepresented as being the original software.
 
* Free for non-commercial use only.  For commercial use, explicit approval 
must be requested by contacting the Authors.

* If you use the code in your work, you must acknowledge it

* Modifications of the source code must also be released as open source

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/


#ifndef _nodetable_h
#define _nodetable_h

#include <vector>
#include <map>

using namespace std;

/** \class NodeTable
 * \brief Class that contains node references.
 *
 * Class that contains references to all of the nodes in the vascular structure
 * as well as the relationship between the nodes (the nodes ancestor and siblings).
 * It also contains information about the flow, radius, and reduced resistance.
 * It essentially describes the generated vascular structure.
 *
 * The NodeTable is essentially a vector of arrays.  The arrays are a series of doubles
 * that correspond to various properties of the nodes in the vascular structure.  
 * The values that are stored include the node type, position, parent node, leftRatio,
 * rightRatio, flow, leftChild, rightChild in that order.
 *
 */
class NodeTable {

public:
	//constants used for the type field of a node
	static size_t TERM;
	static size_t ROOT;
	static size_t BIF;
	static size_t FIELDS;

	//copy of entries in the node table before they were modified (used for undos)
	map<int, double*> originalEntries;
  size_t length;
	bool prepareUndo;

	//current set of nodes
	vector<double *> nodes;
	
	NodeTable();

	//undo utilities
	void startUndo();
	void stopUndo();
	void clearUndo();
	void applyUndo();
	
	
	//setters/getters
  size_t getType(size_t index);
	void setType(size_t index, size_t type);
	double* getPos(size_t index);
	void setPos(size_t index, double* source);
  size_t getParent(size_t index);
	void setParent(size_t index, size_t parent);
	double getLeftRatio(size_t index);
	void setLeftRatio(size_t index, double ratio);
	double getRightRatio(size_t index);
	void setRightRatio(size_t index, double ratio);
	double getFlow(size_t index);
	void setFlow(size_t index, double flow);
  size_t getLeftChild(size_t index);
	void setLeftChild(size_t index, size_t id);
  size_t getRightChild(size_t index);
	void setRightChild(size_t index, size_t id);
	double getRadius(size_t index);
	void setRadius(size_t index, double radius);
	double getReducedResistance(size_t index);
	void setReducedResistence(size_t index, double resistance);
	
	//adds a new node to the node table (makes a copy of *node
	void addNode(double* node);
	//sets the node specified by index to *node
	void setNode(size_t index, double* node);
	//adds a new node to the node table
	void addNode(size_t type, double* pos, size_t parent, double leftRatio, double rightRatio, double flow, size_t leftChild, size_t rightChild);

	//makes a copy of the current node table
	NodeTable copy();
};
#endif
