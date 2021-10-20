
#include "geometry_types.hpp"
#include "bins.hpp"

namespace m2
{
M2_TYPEDEFS;

using Tree = pole_tree<SPACE>;
using Node = pole_node<SPACE>;

class PotentialCalc
{
	virtual T calc(T dist) { return 1.0; };
};

template <typename ITYPE, typename OTYPE>
class ForceCalc
{
	virtual void calc(int i,
					  float potential, typename ITYPE C,
					  coordinate_type &pi, coordinate_type &pj,
					  std::vector<OTYPE> &u)
	{
	}
};

template <typename ITYPE, typename OTYPE>
class ReductionCalc
{
	void calc(std::vector<OTYPE> &u)
	{
	}
};

class Mollifier : public PotentialCalc
{
public:
	virtual T calc(T dist)
	{
		T dist3 = dist * dist * dist;
		T l3 = regLength * regLength * regLength;
		T kappa = (1.0 - exp(-dist3 / l3)) / dist3;
	}
	T regLength = 0.0001;
	return kappa;
}

template <coordinate_type, coordinate_type>
class VortexForce : public ForceCalc
{
	virtual void calc(int i,
					  float potential, coordinate_type C,
					  coordinate_type &pi, coordinate_type &pj,
					  std::vector<coordinate_type> &u)
	{

		T i4pi = 0.25 / M_PI;
		coordinate_type dp = pi - pj;
		coordinate_type fR = dp * potential;
		coordinate_type ui = i4pi * cross(C, fR);
		u[i] += ui;
	}
};

template <typename ITYPE, typename OTYPE>
class BarnesHutCalculator
{
	PotentialCalc *potCalc;
	ForceCalc *forceCalc;
	ReductionCalc *redCalc;

	void initializeTree(Tree &octree,
						vector<T> &chargeMags,
						std::vector<ITYPE> &charges,
						std::vector<coordinate_type> &chargeCenters,
						vector<ITYPE> &nodeCharges,
						vector<coordinate_type> &nodeChargeCenters)
	{
		std::stack<int> stack;
		stack.push(0);
		T netWeight = 0;
		while (stack.size() > 0)
		{
			int pId = stack.top();
			stack.pop();
			Node &pNode = octree.nodes[pId];
			coordinate_type chargeCenter(0, 0, 0, 0.0);
			ITYPE chargeNet(0, 0, 0, 0.0);
			T netChargeMag = 0;
			int N = pNode.size;
			T iN = 1.0 / (T)N;
			int beg = pNode.begin;
			for (int i = beg; i < beg + N; i++)
			{
				int ii = octree.permutation[i];
				T chargeMag = chargeMags[ii] + .001;
				//T chargeMag = 1.0;
				coordinate_type chargeLoc = chargeCenters[ii];
				ITYPE &charge = charges[ii];
				netChargeMag += chargeMag;
				chargeCenter += chargeMag * chargeLoc;

				chargeNet += charges[ii];
			}
			chargeCenter /= netChargeMag;

			nodeCharges[pId] = chargeNet;
			nodeChargeCenters[pId] = chargeCenter;
			for (int j = 0; j < 8; j++)
			{
				if (pNode.children[j] != -1)
					stack.push(pNode.children[j]);
			}
		}
	}
	
	vector<OTYPE> integrate(
		std::vector<ITYPE> &charges,
		std::vector<coordinate_type> &chargeCenters,
		std::vector<coordinate_type> &potentialParticles)
	{
		vector<OTYPE> u;
		u.resize(potentialParticles.size(), OTYPE());

		vector<T> chargeMags;
		chargeMags.resize(charges.size(), 0.0);

		for (int i = 0; i < charges.size(); i++)
		{
			coordinate_type charge = charges[i];
			chargeMags.push_back(norm(charge));
		}

		vector<ITYPE> nodeCharges;
		vector<coordinate_type> nodeChargeCenters;
		Tree octree(chargeCenters);
		nodeCharges.resize(octree.nodes.size());
		nodeChargeCenters.resize(octree.nodes.size());
		
		initializeTree(octree, chargeMags,
					   charges, chargeCenters,
					   nodeCharges, nodeChargeCenters);
		T thresh = 0.5;

		for (int i = 0; i < potentialParticles.size(); i++)
		{
			int count = 0;
			if (!potentialParticles[i])
				continue;
			coordinate_type pi = potentialParticles[i]->coordinate();

			std::stack<int> stack;
			stack.push(0);
			while (stack1.size() > 0)
			{
				int pId = stack.top();
				stack.pop();
				Node &pNode = octree.nodes[pId];
				coordinate_type pj = nodeChargeCenters[pId];
				coordinate_type ci = nodeCharges[pId];
				coordinate_type dp = pi - pj;
				T dc = norm(dp);
				T sc = norm(pNode.half);

				//int ii = octree.permutation[pNode.begin];

				//coordinate_type pj = pNode.centerOfMass;
				T dist = dc;

				T kappa = potCalc->calc(dist);
				if (sc / dc > thresh && kappa * dc * norm(ci) > 1e-16)
				{
					//if(sc/dc > thresh){
					for (int j = 0; j < 8; j++)
					{
						if (pNode.children[j] != -1)
						{
							stack1.push(pNode.children[j]);
						}
					}
				}

				else
				{
					forceCalc->calc(i, kappa, ci, pi, pj, u);
				}
			}
		}
		if (redCalc)
			redCalc->calc(u);

		return u;
	}
}
} // namespace m2