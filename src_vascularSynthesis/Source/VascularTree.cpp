// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com
/*=========================================================================
 
	Program: VascuSynth
	Module: $RCSfile: VascularTree.cpp,v $
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

#include <cmath>
#include <iostream>
#include <random>
#include <time.h>
#include <chrono>
#include <vector>
#include <stdio.h>
#include <algorithm>

#include "VascularTree.h"
#include "OxygenationMap.h"
#include "NodeTable.h"

#define PI 3.1415926535897


using namespace std;

/**  
 * 	Constructor
 */				
VascularTree::VascularTree(OxygenationMap * _oxMap, double* _perf, double _Pperf, double _Pterm, double _Qperf, double _rho, double _gamma, double _lambda, double _mu, double _minDistance, int _numNodes, double _voxelWidth, int _closestNeighbours, double _tempreture){
	oxMap = _oxMap;
	perf = new double[3]; perf[0] = _perf[0]; perf[1] = _perf[1];perf[2] = _perf[2];
	Pperf = _Pperf;
	Pterm = _Pterm;
	Qperf = _Qperf;
	rho = _rho;
	gamma = _gamma;
	lambda = _lambda;
	mu = _mu;
	minDistance = _minDistance;
	mapVoxelWidth = _voxelWidth;
	Qterm = _Qperf/_numNodes;
	numNodes = _numNodes;
	
	closestNeighbours = _closestNeighbours;
	
	tempreture = _tempreture;
	
	nt.addNode(NodeTable::ROOT, perf, -1, 1, 1, Qperf, -1, -1);
}

/** 
 * Calculates the distance between to nodes in the node table.
 */
double VascularTree::distance(size_t from, size_t to){
	double* fromPos = nt.getPos(from);
	double* toPos = nt.getPos(to);
	
	return  sqrt(
				pow(fromPos[0] - toPos[0], 2) + 
				pow(fromPos[1] - toPos[1], 2) +
				pow(fromPos[2] - toPos[2], 2))*mapVoxelWidth;
}

/**
 * Calculates the reduced resistence of segment at id.
 */
void VascularTree::calculateReducedResistence(size_t id){
	if(nt.getType(id) == NodeTable::TERM){
		double acc = (8.0 * rho * distance(id, nt.getParent(id))/PI); 
		nt.setReducedResistence(id, acc);
	} else {
		double acc = 0;
		acc += pow(nt.getLeftRatio(id), 4) / nt.getReducedResistance(nt.getLeftChild(id));
		acc += pow(nt.getRightRatio(id), 4) / nt.getReducedResistance(nt.getRightChild(id));
		acc = 1.0 / acc;
		acc += (8.0 * rho * distance(id, nt.getParent(id))/PI);
		nt.setReducedResistence(id, acc);
	}
}

/**
 * Calculates the ratio of radii of the segment at id.
 */
void VascularTree::calculateRatios(size_t id){
  size_t left = nt.getLeftChild(id);
  size_t right = nt.getRightChild(id);
	
	double left_over_right = (nt.getFlow(left)*nt.getReducedResistance(left)) / (nt.getFlow(right) * nt.getReducedResistance(right));
	left_over_right = pow(left_over_right, 0.25);
	
	nt.setLeftRatio(id, pow((1 + pow(left_over_right, -gamma)), (-1.0)/gamma));
	nt.setRightRatio(id, pow((1 + pow(left_over_right, gamma)), (-1.0)/gamma));
}

/**
 * Updates the tree at the bifurication point at id.
 */
void VascularTree::updateAtBifurication(size_t id, size_t newChild){
	if(nt.getType(id) != NodeTable::ROOT){		
		calculateReducedResistence(newChild);
		calculateRatios(id);
		
		updateAtBifurication(nt.getParent(id), id);
	} else {
		calculateReducedResistence(newChild);
	}
}

/**
 * Calculates the radii throughout the tree.
 */
void VascularTree::calculateRadius(){
	//the child of the root (perferation) node defines the root segment
  size_t rootChild = nt.getLeftChild(0);
	
	double radius = (nt.getFlow(rootChild) * nt.getReducedResistance(rootChild)) / (Pperf - Pterm);
	radius = pow(radius, 0.25);
	nt.setRadius(rootChild, radius);
	
	calculateRadius(rootChild);
}

/**
 * Calculates the radius at id.
 */
void VascularTree::calculateRadius(size_t id){
	if(nt.getType(id) == NodeTable::TERM)
		return;
	
  size_t left = nt.getLeftChild(id);
  size_t right = nt.getRightChild(id);
	
	nt.setRadius(left, nt.getRadius(id)*nt.getLeftRatio(id));
	nt.setRadius(right, nt.getRadius(id)*nt.getRightRatio(id));
	
	calculateRadius(left);
	calculateRadius(right);
}

/**
 * Calculates the fitness.
 */
double VascularTree::calculateFitness(){
	
    size_t size = nt.nodes.size();
    double acc = 0;
    
	for(size_t i = 1; i < size; i++){
		acc += pow(distance(i, nt.getParent(i)), mu) * pow(nt.getRadius(i), lambda);
	}
	
	return acc;
}

/**
 * When used by local optimization, ignored is the segment to connect to
 * otherwise it should be -1;
 */
bool VascularTree::validateCandidate(double* x0, int ignored){
    int size = nt.nodes.size();
    
	for(int i = 1; i < size; i++){
		if(i != ignored){
			double distance = pointSegmentDistance(x0, i);
			if(distance < minDistance)
				return false;
		}
	}
	
	return true;
}

/**
 * Connects the candidate node to a segment through the
 * bifurcation point bifPoint.
 */
void VascularTree::connectPoint(double* point, size_t segment, double* bifPoint){
	if(nt.getType(segment) == NodeTable::ROOT){
		nt.addNode(NodeTable::TERM, point, segment, 1, 1, Qterm, -1, -1);
		nt.setLeftChild(0, 1);
		nt.setRightChild(0, 1);
		calculateReducedResistence(1);
	} else {
		/* *--- <I_seg> --- *
		//			
		//			becomes
		//
		//	*--- <I_bif> --- * ---- < I_con > --- *
		//					 \
		//					  \
		//					  <I_new > 
		//						\
		//						 \
		//						  * point[]
		//
		//
		//  where I_sec is replaced with I_con
		//
		//
		// 	I_con = I_seg with the exception of parent (which is set to I_biff
		//  I_seg's parent updates its child to be I_biff
		//  I_new = is a new segment wich is trivially built
		//  I_biff is built using I_new and I_con */
        
    size_t biffId = nt.nodes.size();
    size_t newId = biffId + 1;
		
    size_t oldParent = nt.getParent(segment);
		
		nt.setParent(segment, biffId);
		if(nt.getLeftChild(oldParent) == segment)
			nt.setLeftChild(oldParent, biffId);
		if(nt.getRightChild(oldParent) == segment)
			nt.setRightChild(oldParent, biffId);
		
		if(oldParent > 0)
			incrementFlow(oldParent, Qterm);
		
		nt.addNode(NodeTable::BIF, bifPoint, oldParent, 1, 1, nt.getFlow(segment) + Qterm, segment, newId);
		nt.addNode(NodeTable::TERM, point, biffId, 1, 1, Qterm, -1, -1);
		
		calculateReducedResistence(segment);
		updateAtBifurication(biffId, newId);
	}
}

/**
 * Updates the flow throughout the tree.
 */
void VascularTree::incrementFlow(size_t parent, double Qterm){
	nt.setFlow(parent, nt.getFlow(parent)+Qterm);
	if(nt.getParent(parent) > 0)
		incrementFlow(nt.getParent(parent), Qterm);
}

/**
 * Returns the distance between a point and a segment.
 */
double VascularTree::pointSegmentDistance(double* x0, size_t segment){
	double *pos1 = nt.getPos(segment);
	double *pos2 = nt.getPos(nt.getParent(segment));

	double t = -((pos2[0] - pos1[0])*(pos1[0] - x0[0]) + 
				 (pos2[1] - pos1[1])*(pos1[1] - x0[1]) + 
				 (pos2[2] - pos1[2])*(pos1[2] - x0[2])) / 
				((pos2[0] - pos1[0])*(pos2[0] - pos1[0]) + 
				 (pos2[1] - pos1[1])*(pos2[1] - pos1[1]) +
				 (pos2[2] - pos1[2])*(pos2[2] - pos1[2]));
	
	if(t < 0 || t > 1){
		double d1 = pow(pos1[0] - x0[0], 2) + pow(pos1[1] - x0[1], 2) + pow(pos1[2] - x0[2], 2);
		double d2 = pow(pos2[0] - x0[0], 2) + pow(pos2[1] - x0[1], 2) + pow(pos2[2] - x0[2], 2);
        
		if(d1 < d2)
			return pow(d1, 0.5);
		else
			return pow(d2, 0.5);
	}
	else{
		return  pow(pow((((pos2[0] - pos1[0])*t + pos1[0]) - x0[0]), 2) + 
						 pow((((pos2[1] - pos1[1])*t + pos1[1]) - x0[1]), 2) + 
						 pow((((pos2[2] - pos1[2])*t + pos1[2]) - x0[2]), 2), 0.5);
	}
}

/**
 * Optimizes the location of a bifurication point for terminal node 
 * point and segment 'segment'
 */
double* VascularTree::localOptimization(double * point, size_t segment, int steps){
	double bestFitness = 1e200;
	//double tempreture = 1e-7;//1e-3;//1e-10;//0.0003;

	double bif[3];
	double perf[3];
	perf[0] = nt.getPos(nt.getParent(segment))[0];
	perf[1] = nt.getPos(nt.getParent(segment))[1];
	perf[2] = nt.getPos(nt.getParent(segment))[2];
	double con[3];
	con[0] = nt.getPos(segment)[0];
	con[1] = nt.getPos(segment)[1];
	con[2] = nt.getPos(segment)[2];
	
	bif[0] = ((con[0] - perf[0])/2.0 + perf[0]);
	bif[1] = ((con[1] - perf[1])/2.0 + perf[1]);
	bif[2] = ((con[2] - perf[2])/2.0 + perf[2]);
	
	double stepSize = (((perf[0]+con[0]+point[0])/3.0) - bif[0]+
					   ((perf[1]+con[1]+point[1])/3.0) - bif[1] +
					   ((perf[2]+con[2]+point[2])/3.0) - bif[2])* 2.0 / steps;
					   
	int sampleNumber = 100;
    vector<vector<double>> sampling(sampleNumber, vector<double>(2));

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0,1.0);
    
    for(int i=0; i<sampleNumber; i++){
    	sampling[i][0] = dis(gen);
    	sampling[i][1] = dis(gen);
    }
	
	//makesure point is visible from bif - otherwise return null
	if(!oxMap->visible(bif, point) || !inVolume(bif))
		return NULL;
	
	nt.startUndo();

	connectPoint(point, segment, bif);
	nt.applyUndo();
		
	double localBest[3];
	double test[3];
	vector<double> fitnessDistribution(sampleNumber,0);
	vector<double> fitnessRecord(sampleNumber,1e100);
	vector<vector<double>> testBif(sampleNumber, vector<double>(3,0));
	
	for(int i = 0; i<sampleNumber; i++){
		double r1 = sampling[i][0];
		double r2 = sampling[i][1];
		
		test[0] = (1 - sqrt(r1)) * perf[0] + (sqrt(r1)*(1 - r2)) * con[0] + (r2*sqrt(r1)) * point[0];
		test[1] = (1 - sqrt(r1)) * perf[1] + (sqrt(r1)*(1 - r2)) * con[1] + (r2*sqrt(r1)) * point[1];
		test[2] = (1 - sqrt(r1)) * perf[2] + (sqrt(r1)*(1 - r2)) * con[2] + (r2*sqrt(r1)) * point[2];
		
		testBif[i][0] = test[0];
		testBif[i][1] = test[1];
		testBif[i][2] = test[2];
		
		if(r1==0 || (r1==1 && r2==0) || (r1==1 && r2==1)) // bif should not be perf, con or point
			continue;
			
		double stepsCur = (test[0] - bif[0] + test[1] - bif[1] + test[2] - bif[2] ) * 1.0 / stepSize;
		if(stepsCur>steps)
			continue;
		
		if(inVolume(test) &&oxMap->visible(perf, test) && oxMap->visible(con, test) && oxMap->visible(point, test) && validateCandidate(test, segment)){
			connectPoint(point, segment, test);
			calculateRadius();
			double fitness = calculateFitness();
			
			//fitnessDistribution[i] = exp(-1.0*fitness / tempreture);
			fitnessRecord[i] = fitness;
			
		}
		nt.applyUndo();
	}
	
	//for(int i=0;i<sampleNumber;i++)
		//cout<<fitnessDistribution[i]<<", ";
	
	default_random_engine generator;
	//discrete_distribution<int> distribution (fitnessDistribution.begin(),fitnessDistribution.end());
	
	int bifIndex;
	//int bifIndexTmp;
	//if(1)//(tempreture <1e-4)
		bifIndex = min_element(fitnessRecord.begin(),fitnessRecord.end()) - fitnessRecord.begin();
	/*else{
		vector<int> count(sampleNumber,0);
		for(int k = 0; k<sampleNumber; k++)
			count[distribution(generator)]++;
		bifIndex = max_element(count.begin(), count.end()) - count.begin();
		//bifIndexTmp = distribution(generator); 
	}*/
	for(int i=0; i< sampleNumber;i++)
		fitnessDistribution[i] = exp(-1.0*(fitnessRecord[i]-fitnessRecord[bifIndex])/tempreture);
	discrete_distribution<int> distribution (fitnessDistribution.begin(),fitnessDistribution.end());
	//bifIndex = distribution(generator);
	localBest[0] = testBif[bifIndex][0]; localBest[1] = testBif[bifIndex][1]; localBest[2] = testBif[bifIndex][2]; 
	//bestFitness = fitnessRecord[bifIndex];
	
	double ratio = 0.5; // control the standard deviation
	if(dis(gen)<=0.49){ // control the mean
		localBest[0] = (1-ratio)*localBest[0] + ratio*perf[0];
		localBest[1] = (1-ratio)*localBest[1] + ratio*perf[1];
		localBest[2] = (1-ratio)*localBest[2] + ratio*perf[2];
		connectPoint(point, segment, localBest);
		calculateRadius();
		bestFitness = calculateFitness();
		nt.applyUndo();
	}
	else{
		localBest[0] = (1-ratio)*localBest[0] + ratio*con[0];
		localBest[1] = (1-ratio)*localBest[1] + ratio*con[1];
		localBest[2] = (1-ratio)*localBest[2] + ratio*con[2];
		connectPoint(point, segment, localBest);
		calculateRadius();
		bestFitness = calculateFitness();
		nt.applyUndo();
	}
		
	//cout<<endl<<min_element(fitnessRecord.begin(),fitnessRecord.end()) - fitnessRecord.begin()<<endl;
	//cout<<bifIndex<<", "<<bifIndexTmp<<endl;
	//cin.get();
    //bestFitness = fitnessRecord[bifIndex];

    bif[0] = localBest[0]; bif[1] = localBest[1]; bif[2] = localBest[2];

	/*for(int i = 0; i < steps; i++){
		localBest[0] = bif[0]; localBest[1] = bif[1]; localBest[2] = bif[2];
		//try the neighbours of bif (6 connected) as possible new bif points
		//		for a bif point to be valid:
		//			- perf must be visible from bif
		//			- con must be visible from bif
		//			- point must be visilbe from bif
		//			- bif must be a valid candidate (using segment as the ignored)
		//
		//	if the current bif is ever a local optima -search is terminated
		
		vector<double> fitnessDistribution(6,1e-10); // Gibbs Distribution for six directions
		vector<double*> testBif(6);
		vector<double> fitnessRecord(6,1e200);
		
		//up
		test[0] = bif[0]+stepSize; test[1] = bif[1]; test[2] = bif[2];
		double test0[3];
		test0[0] = test[0]; test0[1] = test[1]; test0[2] = test[2];
		testBif[0] = test0;
		
		if(inVolume(test) &&oxMap->visible(perf, test) && oxMap->visible(con, test) && oxMap->visible(point, test) && validateCandidate(test, segment)){
			connectPoint(point, segment, test);
			calculateRadius();
			double fitness = calculateFitness();
			
			fitnessDistribution[0] = exp(-1*fitness / tempreture);
			fitnessRecord[0] = fitness;
			//if(fitness < bestFitness){
			//	localBest[0] = test[0]; localBest[1] = test[1]; localBest[2] = test[2];
			//	bestFitness = fitness;
			//}
		}
		nt.applyUndo();
		
		//down
		test[0] = bif[0]-stepSize; test[1] = bif[1]; test[2] = bif[2];
		double test1[3];
		test1[0] = test[0]; test1[1] = test[1]; test1[2] = test[2];
		testBif[1] = test1;
		
		if(inVolume(test) &&oxMap->visible(perf, test) && oxMap->visible(con, test) && oxMap->visible(point, test) && validateCandidate(test, segment)){
			connectPoint(point, segment, test);
			calculateRadius();
			double fitness = calculateFitness();
			
			fitnessDistribution[1] = exp(-1*fitness / tempreture);
			fitnessRecord[1] = fitness;
			//if(fitness < bestFitness){
			//	localBest[0] = test[0]; localBest[1] = test[1]; localBest[2] = test[2];
			//	bestFitness = fitness;
			//}
		}
		nt.applyUndo();

		//left
		test[0] = bif[0]; test[1] = bif[1]+stepSize; test[2] = bif[2];
		double test2[3];
		test2[0] = test[0]; test2[1] = test[1]; test2[3] = test[3];
		testBif[2] = test2;
		
		if(inVolume(test) &&oxMap->visible(perf, test) && oxMap->visible(con, test) && oxMap->visible(point, test) && validateCandidate(test, segment)){
			connectPoint(point, segment, test);
			calculateRadius();
			double fitness = calculateFitness();
			
			fitnessDistribution[2] = exp(-1*fitness / tempreture);
			fitnessRecord[2] = fitness;
			//if(fitness < bestFitness){
			//	localBest[0] = test[0]; localBest[1] = test[1]; localBest[2] = test[2];
			//	bestFitness = fitness;
			//}
		}
		nt.applyUndo();

		//right
		test[0] = bif[0]; test[1] = bif[1]-stepSize; test[2] = bif[2];
		double test3[3];
		test3[0] = test[0]; test3[1] = test[1]; test3[2] = test[2];
		testBif[3] = test3;
		
		if(inVolume(test) &&oxMap->visible(perf, test) && oxMap->visible(con, test) && oxMap->visible(point, test) && validateCandidate(test, segment)){
			connectPoint(point, segment, test);
			calculateRadius();
			double fitness = calculateFitness();
			
			fitnessDistribution[3] = exp(-1*fitness / tempreture);
			fitnessRecord[3] = fitness;
			//if(fitness < bestFitness){
			//	localBest[0] = test[0]; localBest[1] = test[1]; localBest[2] = test[2];
			//	bestFitness = fitness;
			//}
		}
		nt.applyUndo();

		//forward
		test[0] = bif[0]; test[1] = bif[1]; test[2] = bif[2]+stepSize;
		double test4[3];
		test4[0] = test[0]; test4[1] = test[1]; test4[2] = test[2];
		testBif[4] = test4;
		
		if(inVolume(test) &&oxMap->visible(perf, test) && oxMap->visible(con, test) && oxMap->visible(point, test) && validateCandidate(test, segment)){
			connectPoint(point, segment, test);
			calculateRadius();
			double fitness = calculateFitness();
			
			fitnessDistribution[4] = exp(-1*fitness / tempreture);
			fitnessRecord[4] = fitness;
			//if(fitness < bestFitness){
			//	localBest[0] = test[0]; localBest[1] = test[1]; localBest[2] = test[2];
			//	bestFitness = fitness;
			//}
		}
		nt.applyUndo();

		//back
		test[0] = bif[0]; test[1] = bif[1]; test[2] = bif[2]-stepSize;
		double test5[3];
		test5[0] = test[0]; test5[1] = test[1]; test5[2] = test[2];
		testBif[5] = test5;
		
		if(inVolume(test) &&oxMap->visible(perf, test) && oxMap->visible(con, test) && oxMap->visible(point, test) && validateCandidate(test, segment)){
			connectPoint(point, segment, test);
			calculateRadius();
			double fitness = calculateFitness();
			
			fitnessDistribution[5] = exp(-1*fitness / tempreture);
			fitnessRecord[5] = fitness;
			//if(fitness < bestFitness){
			//	localBest[0] = test[0]; localBest[1] = test[1]; localBest[2] = test[2];
			//	bestFitness = fitness;
			//}
		}
		nt.applyUndo();
		
		//if(localBest[0] != bif[0] || localBest[1] != bif[1] || localBest[2] != bif[2]){
		//	bif[0] = localBest[0]; bif[1] = localBest[1]; bif[2] = localBest[2];
		//} else {
		//	break;
		//}
		
		//use this distribution to sample
		default_random_engine generator;
		discrete_distribution<int> distribution (fitnessDistribution.begin(),fitnessDistribution.end());
		
		int direction = distribution(generator);
		localBest[0] = testBif[direction][0]; localBest[1] = testBif[direction][1]; localBest[2] = testBif[direction][2]; // here "localBest" is not the local minima anymore
        bestFitness = fitnessRecord[direction]; //also "bestFitness" is not the best anymore
        
		bif[0] = localBest[0]; bif[1] = localBest[1]; bif[2] = localBest[2];
	}*/
	nt.clearUndo();
	nt.stopUndo();
	
	double *ret = new double[4];
	ret[0] = bif[0]; ret[1] = bif[1]; ret[2] = bif[2];
	ret[3] = bestFitness;
	return ret;
}

/**
 * Determines if point is in the volume.
 */
bool VascularTree::inVolume(double* point){
	if(point[0] < 0 || point[0] >= oxMap->dim[0])
		return false;
	if(point[1] < 0 || point[1] >= oxMap->dim[1])
		return false;
	if(point[2] < 0 || point[2] >= oxMap->dim[2])
		return false;
	
	return true;
}

/**
 *  Connects the candidate node to the closestNeighbour segments
 *  with an optimized bifurcation location (originally the midpoint of the segment.
 *  The conceptually connected segment that yields the smallest objective
 *  function is elected.
 */
bool VascularTree::connectCandidate(double* point, int steps){
	
	//cout << "connect candidate" << endl;
	
	if(!validateCandidate(point, -1)) {
		return false;	//candiate is too close to an existing segment
	}

	if(nt.nodes.size() == 1){
		
		if(!oxMap->visible(nt.getPos(0), point)) {
			return false;			
		}
		
		connectPoint(point, 0, NULL);
		return true;
	}
	
	double best[3];
	bool foundSolution = false;
	double bestFitness = 1e200;
  size_t bestSegment = 0;
	
	map<int, double> distances;
  size_t i;
  size_t j;
	
	//determine the distance between the point and each segment
  size_t size = nt.nodes.size();
    for(i = 1; i < size; i++){
		double d = pointSegmentDistance(point, i);
		distances[i] = d;
	}
	
  size_t numCandidates = distances.size();
	size_t *candidateSegments = new size_t[numCandidates];
	
	//sort the segments by distance
	for(j = 0; j < numCandidates; j++){
		int closest = -1;
		double distance = 1e200;
		
		map<int, double>::iterator itr;
		for(itr = distances.begin(); itr != distances.end(); ++itr){
			double d = itr->second;
			if(d < distance){
				distance = d;
				closest = itr->first;
			}
		}
		if(closest >= 1)
			distances.erase(closest);
    if (closest == -1)
    {
      delete[] candidateSegments;
      return false;
    }
		candidateSegments[j] = closest;
	}
	
	//try the first 'closestNeighbours' segments
  size_t count = 0;
	double* test;
	for(j = 0; j < numCandidates && count < closestNeighbours; j++){
		test = localOptimization(point, candidateSegments[j], steps);
		if(test != NULL){
			count++;
			if(test[3] < bestFitness){
				best[0] = test[0]; best[1] = test[1]; best[2] = test[2];
				bestSegment = candidateSegments[j];
				bestFitness = test[3];
				foundSolution = true;
			}
		}
		delete test;
	}
	
	delete[] candidateSegments;
    

	if(!foundSolution)
		return false; //could not connect candidate to ANY segment
	
	connectPoint(point, bestSegment, best);
	
	return true;
}

/**
 * iteratively builds a tree by selecting candidate nodes based on an oxygenation demand map
 * and then connecting that candidate node to a series of closest segments at the midpoint of the
 * segment.  The resulting new segment that yields the smallest objective function is selected and 
 * added to the tree.  The oxygenation map is updated - reducing values proximal to the candidate node.
 *
 * finally the radii are calcualted recursively by calculating the ratio of the radius of the child over parent,
 * and multiplying that by the radius of the parent to ascertain the radius.
 */
void VascularTree::buildTree(){
	
	int count = 0;
	int count_cur;
	double cand[3];
    size_t term[3];
	while(count < numNodes){
		
		count_cur = count;
		double sum = oxMap->sum();
		oxMap->candidate(sum, term);
		cand[0] = term[0]; cand[1] = term[1]; cand[2] = term[2];
		if(connectCandidate(cand, 20)){
			count++;
			oxMap->applyCandidate(term);
		}
		if(count_cur != count){
			time_t my_time = time(NULL); 
			cout<<"system time: "<<ctime(&my_time)<<"count: " << count << endl;
		}
	}
	
	calculateRadius();
	
}
