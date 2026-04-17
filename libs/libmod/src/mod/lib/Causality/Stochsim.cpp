#include "Stochsim.hpp"

#include <mod/lib/Random.hpp>

#include <jla_boost/graph/PairToRangeAdaptor.hpp>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <algorithm>
#include <iostream>
#include <cmath>

namespace mod::lib::Causality {
namespace {

using PropensityEntry = std::pair<int, double>;

double reactionPropensity(lib::DG::HyperVertex e, const Marking &m) {
	const petri::Transition t = m.getNet().getTransition(e);
	const auto &marking = m.getMarking();
	assert(marking.isEnabled(t));
	const auto &net = m.getNet().getNet();
	const auto &g = net.getGraph();
	const auto vt = net.vertexFromTransition(t);
	double res = 1.0;
	for(const auto eIn: asRange(in_edges(vt, g))) {
		const auto vIn = source(eIn, g);
		assert(g[vIn].kind == petri::Net::Kind::Place);
		const int c = marking[net.placeFromVertex(vIn)];
		const int w = g[eIn];
		switch(w) {
		case 1:
			res *= c;
			break;
		case 2:
			res *= c * (c - 1) / 2;
			break;
		default:
			res *= boost::math::binomial_coefficient<double>(c, w);
			break;
		}
	}
	return res;
}

std::vector<PropensityEntry> computePropensities(
		const lib::DG::Hyper &dg,
		const Marking &m,
		const std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> &inputRate,
		const std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> &reactionRate,
		const std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> &outputRate,
		std::vector<double> &cachedInputRates,
		std::vector<double> &cachedRates) {
	const auto &dgGraph = dg.getGraph();

	std::vector<PropensityEntry> propensities; // .first: non-negative==reaction/output, negative: -input - 1
	propensities.reserve(num_vertices(dgGraph));

	for(const auto e: m.getAllEnabled()) {
		const auto idx = get(boost::vertex_index_t(), dgGraph, e);
		assert(idx < cachedRates.size());
		double r = cachedRates[idx];
		if(r < 0) {
			if(reactionRate) {
				bool cache;
				std::tie(r, cache) = reactionRate(dg, e);
				assert(r >= 0);
				if(cache) cachedRates[idx] = r;
			} else {
				cachedRates[idx] = r = 1.0;
			}
		}
		if(r != 0) propensities.emplace_back(idx, r * reactionPropensity(e, m));
	}
	for(const auto v: m.getNonZeroPlaces()) {
		const auto idx = get(boost::vertex_index_t(), dgGraph, v);
		assert(idx < cachedRates.size());
		double r = cachedRates[idx];
		if(r < 0) {
			if(outputRate) {
				bool cache;
				std::tie(r, cache) = outputRate(dg, v);
				assert(r >= 0);
				if(cache) cachedRates[idx] = r;
			} else {
				cachedRates[idx] = r = 0.0;
			}
		}
		if(r != 0) propensities.emplace_back(idx, r * m.getMarking()[m.getNet().getPlace(v)]);
	}
	for(const auto v: asRange(vertices(dgGraph))) {
		if(dgGraph[v].kind != lib::DG::HyperVertexKind::Vertex) continue;
		const auto idx = get(boost::vertex_index_t(), dgGraph, v);
		assert(idx < cachedInputRates.size());
		double r = cachedInputRates[idx];
		if(r < 0) {
			if(inputRate) {
				bool cache;
				std::tie(r, cache) = inputRate(dg, v);
				assert(r >= 0);
				if(cache) cachedInputRates[idx] = r;
			} else {
				cachedInputRates[idx] = r = 0.0;
			}
		}
		if(r != 0) propensities.emplace_back(-idx - 1, r);
	}

	return propensities;
}

std::vector<std::pair<petri::Place,int>> consumed(const lib::DG::Hyper &dg, const Marking &m, int idx) {
    const auto &dgGraph = dg.getGraph();
    if(idx >= 0) {
        const auto v = vertices(dgGraph).first[idx];
        if(dgGraph[v].kind == lib::DG::HyperVertexKind::Edge) {
            const auto &net = m.getNet().getNet();
            const auto t = m.getNet().getTransition(v);
            std::vector<std::pair<petri::Place,int>> result;
            for(const auto &[place, w] : net.consumed(t)) {
                std::pair<petri::Place,int> pair(place,w);
                result.emplace_back(pair);
            }
            return result;
        } else if(dgGraph[v].kind == lib::DG::HyperVertexKind::Vertex) {
            const auto place = m.getNet().getPlace(v);
            return std::vector<std::pair<petri::Place,int>>{{place,1}};
        } else return {};
    } else return {};
}

std::vector<std::pair<petri::Place,int>> produced(const lib::DG::Hyper &dg, const Marking &m, int idx) {
    const auto &dgGraph = dg.getGraph();
    if(idx >= 0) {
        const auto v = vertices(dgGraph).first[idx];
        if(dgGraph[v].kind == lib::DG::HyperVertexKind::Edge) {
            const auto &net = m.getNet().getNet();
            const auto t = m.getNet().getTransition(v);
            std::vector<std::pair<petri::Place,int>> result;
            for(const auto &[place, w] : net.produced(t)) {
                std::pair<petri::Place,int> pair(place,w);
                result.emplace_back(pair);
            }
            return result;
        } else return {};
    } else {
        const auto v = vertices(dgGraph).first[-idx-1];
        const auto place = m.getNet().getPlace(v);
        return std::vector<std::pair<petri::Place,int>>{{place,1}};
    }
}

} // namespace

DrawMassActionFunction::DrawMassActionFunction(
		const lib::DG::Hyper &dg,
		std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> inputRate,
		std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> reactionRate,
		std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> outputRate)
	: dg(dg), inputRate(inputRate), reactionRate(reactionRate), outputRate(outputRate) {
	syncSize();
}

void DrawMassActionFunction::syncSize() {
	const auto &g = dg.getGraph();
	const auto n = num_vertices(g);
	cachedInputRates.resize(n, -1.0);
	cachedRates.resize(n, -1.0);
}

std::tuple<Action, double, bool> DrawMassActionFunction::draw(const Marking &m) {
	return draw_v0(m);
}

std::tuple<Action, double, bool> DrawMassActionFunction::draw_v0(const Marking &m) {
	constexpr bool VERBOSE = true;

	if(VERBOSE) std::cout << __func__ << ":" << __LINE__ << ":" << std::endl;

	const auto &dgGraph = dg.getGraph();
	auto propensities = computePropensities(dg, m, inputRate, reactionRate, outputRate,
	                                        cachedInputRates, cachedRates);

	if(propensities.empty()) {
		if(VERBOSE) std::cout << __func__ << ":" << __LINE__ << ": no actions" << std::endl;
		return {{}, 0.0, false};
	}

	std::vector<double> accPropensities(propensities.size());
	{
		// TODO: when GCC 8 can be dropped, change to the commented code
		//		std::inclusive_scan(propensities.begin(), propensities.end(), accPropensities.begin(),
		//		                    [](double a, std::pair<int, double> b) {
		//			                    return a + b.second;
		//		                    }, 0.0);
		double sum = 0;
		auto out = accPropensities.begin();
		for(const auto &p: propensities) {
			sum += p.second;
			*out = sum;
			++out;
		}
	}

	if(VERBOSE) {
		std::cout << __func__ << ":" << __LINE__ << ": propensities and accPropensities:" << std::endl;
		for(int i = 0; i != propensities.size(); ++i) {
			std::cout << __func__ << ":" << __LINE__
					<< ": " << propensities[i].first
					<< " " << propensities[i].second
					<< " " << accPropensities[i] << std::endl;
		}
	}

	std::uniform_real_distribution<> dist(0, accPropensities.back());
	auto &rng = mod::lib::getRng();
	const double rnd = dist(rng);
	const auto pos = std::lower_bound(accPropensities.begin(), accPropensities.end(), rnd);
	const auto i = pos - accPropensities.begin();
	if(VERBOSE) std::cout << __func__ << ":" << __LINE__ << ": rnd=" << rnd << " i=" << i << std::endl;
	auto actId = propensities[i].first;
	if(actId >= 0) {
		const auto v = vertices(dgGraph).first[actId];
		assert(get(boost::vertex_index_t(), dgGraph, v) == actId);
		if(dgGraph[v].kind == lib::DG::HyperVertexKind::Edge) {
			return {EdgeAction{v}, accPropensities.back(), false};
		} else {
			return {OutputAction{v}, accPropensities.back(), false};
		}
	} else {
		actId = -actId - 1;
		const auto v = vertices(dgGraph).first[actId];
		assert(get(boost::vertex_index_t(), dgGraph, v) == actId);
		assert(dgGraph[v].kind == lib::DG::HyperVertexKind::Vertex);
		return {InputAction{v}, accPropensities.back(), false};
	}
}

// ==============================================================================================

DrawMassActionTauLeapingFunction::DrawMassActionTauLeapingFunction(
		const lib::DG::Hyper &dg,
		std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> inputRate,
		std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> reactionRate,
		std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> outputRate,
		int dc,
		double epsilon)
	: dg(dg), inputRate(inputRate), reactionRate(reactionRate), outputRate(outputRate), dc(dc), epsilon(epsilon) {
	syncSize();
}

void DrawMassActionTauLeapingFunction::syncSize() {
	const auto &g = dg.getGraph();
	const auto n = num_vertices(g);
	cachedInputRates.resize(n, -1.0);
	cachedRates.resize(n, -1.0);
}

std::tuple<Action, double, bool> DrawMassActionTauLeapingFunction::draw(const Marking &m) {
	return draw_v0(m);
}

std::tuple<Action, double, bool> DrawMassActionTauLeapingFunction::draw_v0(const Marking &m) {
	constexpr bool VERBOSE = true;

	if(VERBOSE) std::cout << __func__ << ":" << __LINE__ << ":" << std::endl;

	const auto &dgGraph = dg.getGraph();
	auto propensities = computePropensities(dg, m, inputRate, reactionRate, outputRate,
	                                        cachedInputRates, cachedRates);

	if(propensities.empty()) {
		if(VERBOSE) std::cout << __func__ << ":" << __LINE__ << ": no actions" << std::endl;
		return {{}, 0.0, false};
	}

	if(VERBOSE) {
		std::cout << __func__ << ":" << __LINE__ << ": propensities:" << std::endl;
		for(int i = 0; i != propensities.size(); ++i) {
			std::cout << __func__ << ":" << __LINE__
					<< ": " << propensities[i].first
					<< " " << propensities[i].second << std::endl;
		}
	}

	// idx of the non-critical reactions
	std::vector<int> nonCriticalReactions;
	// propensities for the non-critical reactions
	boost::numeric::ublas::vector<double> nonCriticalPropensities;
	// idx of the critical reactions
    std::vector<int> criticalReactions;
    // propensities for the critical reactions
    std::vector<double> criticalPropensities;
    // temporary buffer for the stoichiometric matrix
    std::vector<boost::numeric::ublas::vector<int>> nonCriticalStoichiometries;

    for(const auto &[idx, propensity] : propensities) {
        // check if the current reactions is critical i.e. if it fully consumes a reactant in less than dc firings
        bool critical = false;
        if(propensity > 0) {
            for(const auto &[place, w] : consumed(dg, m, idx)) {
                if(m.getMarking()[place] / w < dc) {
                    critical = true;
                    break;
                }
            }
	    }

	    if(critical) {
	        criticalReactions.emplace_back(idx);
	        criticalPropensities.emplace_back(propensity);
	    } else {
	        // compute the stoichiometry for the current reaction
	        boost::numeric::ublas::vector<int> deltas(m.getNet().getNet().numPlaces(), 0);
	        for(const auto &[place, w] : consumed(dg, m, idx))
                deltas(place.getId()) -= w;
            for(const auto &[place, w] : produced(dg, m, idx))
                deltas(place.getId()) += w;

            nonCriticalReactions.emplace_back(idx);
            nonCriticalStoichiometries.emplace_back(std::move(deltas));
            nonCriticalPropensities.resize(nonCriticalReactions.size(), true);
            nonCriticalPropensities(nonCriticalReactions.size() - 1) = propensity;
	    }
    }

    // stoichiometric matrix for the non-critical reactions
	boost::numeric::ublas::mapped_matrix<double> stoichiometric(m.getNet().getNet().numPlaces(), nonCriticalReactions.size());

    for(std::size_t reaction = 0; reaction < nonCriticalStoichiometries.size(); ++reaction) {
        const auto &reactionDeltas = nonCriticalStoichiometries[reaction];
        for(unsigned i = 0; i < reactionDeltas.size(); ++i)
            if(reactionDeltas(i) != 0)
                stoichiometric(i, reaction) = static_cast<double>(reactionDeltas(i));
    }

    if(VERBOSE) {
		std::cout << __func__ << ":" << __LINE__ << ": critical reactions:";
		for(int i = 0; i != criticalReactions.size(); ++i)
			std::cout << " " << criticalReactions[i];
	    std::cout << std::endl;
	}

    // compute the highest multiplicity of a reactant of a non-critical reaction for each species
	std::vector<double> gs(stoichiometric.size1());
	for(auto it1 = stoichiometric.begin1(); it1 != stoichiometric.end1(); it1++) {
	    const std::size_t row = it1.index1();
	    for(auto it2 = it1.begin(); it2 != it1.end(); it2++) {
	        const double value = *it2;
	        if(-value > gs[row])
	            gs[row] = -value;
	    }
	}

    // compute the means and variances for the expected firings of any non-critical reaction
	auto sampleMeans = boost::numeric::ublas::prod(stoichiometric, nonCriticalPropensities);
	auto sampleVariances = boost::numeric::ublas::prod(boost::numeric::ublas::element_prod(stoichiometric, stoichiometric), nonCriticalPropensities);

    // compute the maximum tau for which the leap condition holds on the non-critical reactions
    double tau = -1;
    for(const auto v: m.getNonZeroPlaces()) {
        const int amount = m.getMarking()[m.getNet().getPlace(v)];
        const double g = std::max(1.0, gs[m.getNet().getPlace(v).getId()]);

        double sampleMean = sampleMeans(m.getNet().getPlace(v).getId());
        double sampleVariance = sampleVariances(m.getNet().getPlace(v).getId());

        const double newTau = std::min(std::max(amount*epsilon/g,1.0) / std::abs(sampleMean),
                                       std::pow(std::max(amount*epsilon/g,1.0), 2) / sampleVariance);

        if(tau == -1 || newTau < tau) tau = newTau;
    }

    if(VERBOSE) {
        std::cout << __func__ << ":" << __LINE__ << ": tau for non-critical reactions = " << tau << std::endl;
    }

    // compute time till the next criticalReaction
    std::vector<double> accPropensities(criticalPropensities.size());
	{
		// TODO: when GCC 8 can be dropped, change to the commented code
		//		std::inclusive_scan(propensities.begin(), propensities.end(), accPropensities.begin(),
		//		                    [](double a, std::pair<int, double> b) {
		//			                    return a + b;
		//		                    }, 0.0);
		double sum = 0;
		auto out = accPropensities.begin();
		for(const auto &p: criticalPropensities) {
			sum += p;
			*out = sum;
			++out;
		}
	}
	std::uniform_real_distribution<> dist(0,1);
	auto &rng = mod::lib::getRng();
	double timeTillCriticalReaction;
	if(!accPropensities.empty())
	    timeTillCriticalReaction = 1 / accPropensities.back() * std::log(1 / dist(rng));
    else
        timeTillCriticalReaction = -1;

    if(VERBOSE) {
        std::cout << __func__ << ":" << __LINE__ << ": time till next critical reactions = " << timeTillCriticalReaction << std::endl;
    }

    // fire the next critical reaction if it happens withing the tau-interval
    boost::numeric::ublas::vector<double> deltas(m.getNet().getNet().numPlaces(), 0.0);
	if(timeTillCriticalReaction != -1 && (tau == -1 || timeTillCriticalReaction < tau)) {
        tau = timeTillCriticalReaction;

        std::uniform_real_distribution<> dist(0, accPropensities.back());
        auto &rng = mod::lib::getRng();
        const double rnd = dist(rng);
        const auto pos = std::lower_bound(accPropensities.begin(), accPropensities.end(), rnd);
        const auto i = pos - accPropensities.begin();
        if(VERBOSE) {
            std::cout << __func__ << ":" << __LINE__ << ": rnd=" << rnd << " i=" << i << std::endl;
            std::cout << __func__ << ":" << __LINE__ << ": critical reaction fired: " << criticalReactions[i] << std::endl;
        }
        auto actId = criticalReactions[i];

	    for(const auto &[place, w] : consumed(dg, m, actId))
            deltas(place.getId()) -= static_cast<double>(w);
        for(const auto &[place, w] : produced(dg, m, actId))
            deltas(place.getId()) += static_cast<double>(w);
	}

    // compute the expected number of firings for any non-critical reaction
    boost::numeric::ublas::vector<double> firings(nonCriticalPropensities.size());
    for(unsigned reaction = 0; reaction < nonCriticalPropensities.size(); reaction++) {
        std::poisson_distribution<> dist(nonCriticalPropensities(reaction) * tau);
        auto &rng = mod::lib::getRng();
	    firings(reaction) = dist(rng);
    }
    if(VERBOSE) {
        std::cout << __func__ << ":" << __LINE__ << ": non-critical reactions firings:" << std::endl;
        for(int i = 0; i != nonCriticalReactions.size(); ++i)
			std::cout << __func__ << ":" << __LINE__ << ": " << nonCriticalReactions[i] << " " << firings(i) << std::endl;
    }
    boost::numeric::ublas::vector<double> nonCriticalDeltas = boost::numeric::ublas::prod(stoichiometric, firings);
    deltas += nonCriticalDeltas;
    if(VERBOSE) {
        std::cout << __func__ << ":" << __LINE__ << ": deltas:";
        for(int i = 0; i != m.getNet().getNet().numPlaces(); ++i)
			std::cout << " " << deltas(i);
	    std::cout << std::endl;
    }

    // construct an UpdateAction from the deltas
    std::vector<std::pair<lib::DG::HyperVertex, int>> updates;
    for(const auto v: asRange(vertices(dgGraph))) {
        if(dgGraph[v].kind == lib::DG::HyperVertexKind::Vertex) {
            const auto place = m.getNet().getPlace(v);
            const int delta = std::lround(deltas(place.getId()));
            if(delta != 0) updates.emplace_back(v, delta);
        }
    }
    Action action = UpdateAction{std::move(updates)};
    return {action, tau, true};
}

// ==============================================================================================

void Simulator::doIteration() {
	++iteration;
}

} // namespace mod::lib::Causality
