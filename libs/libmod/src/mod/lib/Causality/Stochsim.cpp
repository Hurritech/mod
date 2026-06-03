#include "Stochsim.hpp"

#include <mod/lib/Random.hpp>

#include <jla_boost/graph/PairToRangeAdaptor.hpp>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include <algorithm>
#include <iostream>
#include <cmath>
#include <ctime>

#define USE_TIMER

#ifdef USE_TIMER
    #define INIT_TIMER\
        std::chrono::steady_clock::time_point c_start = std::chrono::steady_clock::now();\
        std::chrono::steady_clock::time_point c_end;\
        int counter = 0;
    #define TIMER do {\
            c_end = std::chrono::steady_clock::now();\
            timings[counter++] += std::chrono::duration_cast<std::chrono::nanoseconds>(c_end - c_start).count();\
            c_start = std::chrono::steady_clock::now();\
        } while(0);
    #define SKIP do { c_start = std::chrono::steady_clock::now(); } while(0);
    #define PRINT_TIMINGS do {\
            if(iteration++ % 1000 == 0) {\
                std::cout << __func__ << __LINE__ << ":";\
                long long sum = 0;\
                for(int i = 0; i < 12; i++)\
                    sum += timings[i];\
                for(int i = 0; i < 12; i++)\
                    std::cout << " " << static_cast<double>(timings[i]) / static_cast<double>(sum);\
                std::cout << std::endl;\
            }\
        } while(0);
#else
    #define INIT_TIMER do { } while(0);
    #define TIMER do { } while(0);
    #define SKIP do { } while(0);
    #define PRINT_TIMINGS do { } while(0);
#endif

namespace mod::lib::Causality {
namespace {

#ifdef USE_TIMER
    static long long timings[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    static int iteration = 0;
#endif

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

double realReactionPropensity(
		lib::DG::HyperVertex e,
		const Marking &m,
		const boost::numeric::ublas::vector<double> &amounts,
		const boost::numeric::ublas::vector<double> *deltas = nullptr) {
	const auto &net = m.getNet().getNet();
	const auto &g = net.getGraph();
	const auto t = m.getNet().getTransition(e);
	const auto vt = net.vertexFromTransition(t);
	double res = 1.0;
	for(const auto eIn: asRange(in_edges(vt, g))) {
		const auto vIn = source(eIn, g);
		assert(g[vIn].kind == petri::Net::Kind::Place);
		const auto place = net.placeFromVertex(vIn);
		double c = amounts(place.getId());
		if(deltas) c += (*deltas)(place.getId());
		c = std::max(0.0, c);
		const int w = g[eIn];
		for(int i = 0; i < w; ++i) {
			if(c <= static_cast<double>(i)) return 0.0;
			res *= (c - static_cast<double>(i)) / static_cast<double>(i + 1);
		}
	}
	return res;
}

bool hasPositiveEntry(const boost::numeric::ublas::vector<double> &v) {
	for(std::size_t i = 0; i < v.size(); ++i)
		if(v(i) > 0.0) return true;
	return false;
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
	constexpr bool VERBOSE = false;

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
		#if defined(__GNUC__) && !defined(__clang__)
            #if __GNUC__ > 8
				std::inclusive_scan(propensities.begin(), propensities.end(), accPropensities.begin(),
				                    [](double a, std::pair<int, double> b) {
					                    return a + b.second;
				                    }, 0.0);
            #else
                // TODO: when GCC 8 can be dropped, remove this
                double sum = 0;
                auto out = accPropensities.begin();
                for(const auto &p: criticalPropensities) {
                    sum += p;
                    *out = sum;
                    ++out;
                }
            #endif
        #endif
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
    INIT_TIMER

	constexpr bool VERBOSE = false;

	if(VERBOSE) std::cout << __func__ << ":" << __LINE__ << ":" << std::endl;

	const auto &dgGraph = dg.getGraph();
	auto propensities = computePropensities(dg, m, inputRate, reactionRate, outputRate,
	                                        cachedInputRates, cachedRates);

    TIMER

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

	SKIP

	// idx of the non-critical reactions
	// std::vector<int> nonCriticalReactions;
	nonCriticalReactions.clear();
	// propensities for the non-critical reactions
	boost::numeric::ublas::vector<double> nonCriticalPropensities;
	// idx of the critical reactions
    // std::vector<int> criticalReactions;
    criticalReactions.clear();
    // propensities for the critical reactions
    // std::vector<double> criticalPropensities;
    criticalPropensities.clear();
    // temporary buffer for the stoichiometric matrix
    std::vector<boost::numeric::ublas::vector<int>> nonCriticalStoichiometries;

    TIMER

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

    TIMER

    // stoichiometric matrix for the non-critical reactions
	boost::numeric::ublas::mapped_matrix<double> stoichiometric(m.getNet().getNet().numPlaces(), nonCriticalReactions.size());

    for(std::size_t reaction = 0; reaction < nonCriticalStoichiometries.size(); ++reaction) {
        const auto &reactionDeltas = nonCriticalStoichiometries[reaction];
        for(unsigned i = 0; i < reactionDeltas.size(); ++i)
            if(reactionDeltas(i) != 0)
                stoichiometric(i, reaction) = static_cast<double>(reactionDeltas(i));
    }

    TIMER

    if(VERBOSE) {
		std::cout << __func__ << ":" << __LINE__ << ": critical reactions:";
		for(int i = 0; i != criticalReactions.size(); ++i)
			std::cout << " " << criticalReactions[i];
	    std::cout << std::endl;
	}

	SKIP

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

	TIMER

    // compute the means and variances for the expected firings of any non-critical reaction
	auto sampleMeans = boost::numeric::ublas::prod(stoichiometric, nonCriticalPropensities);
	auto sampleVariances = boost::numeric::ublas::prod(boost::numeric::ublas::element_prod(stoichiometric, stoichiometric), nonCriticalPropensities);

    TIMER

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

    TIMER

    if(VERBOSE) {
        std::cout << __func__ << ":" << __LINE__ << ": tau for non-critical reactions = " << tau << std::endl;
    }

    SKIP

    // compute time till the next criticalReaction
    // std::vector<double> accPropensities(criticalPropensities.size());
    accPropensities.resize(criticalPropensities.size());
	{
	    #if defined(__GNUC__) && !defined(__clang__)
            #if __GNUC__ > 8
				std::inclusive_scan(criticalPropensities.begin(), criticalPropensities.end(), accPropensities.begin(),
				                    [](double a, double b) {
					                    return a + b;
				                    }, 0.0);
            #else
                // TODO: when GCC 8 can be dropped, remove this
                double sum = 0;
                auto out = accPropensities.begin();
                for(const auto &p: criticalPropensities) {
                    sum += p;
                    *out = sum;
                    ++out;
                }
            #endif
        #endif
	}
	std::uniform_real_distribution<> dist(0,1);
	auto &rng = mod::lib::getRng();
	double timeTillCriticalReaction;
	if(!accPropensities.empty())
	    timeTillCriticalReaction = 1 / accPropensities.back() * std::log(1 / dist(rng));
    else
        timeTillCriticalReaction = -1;

    TIMER

	if(VERBOSE) {
		std::cout << __func__ << ":" << __LINE__ << ": time till next critical reactions = " << timeTillCriticalReaction << std::endl;
	}

	SKIP

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

	TIMER

    // compute the expected number of firings for any non-critical reaction
    boost::numeric::ublas::vector<double> firings(nonCriticalPropensities.size());
    for(unsigned reaction = 0; reaction < nonCriticalPropensities.size(); reaction++) {
        std::poisson_distribution<> dist(nonCriticalPropensities(reaction) * tau);
        auto &rng = mod::lib::getRng();
    	firings(reaction) = dist(rng);
    }

    TIMER

    if(VERBOSE) {
        std::cout << __func__ << ":" << __LINE__ << ": non-critical reactions firings:" << std::endl;
        for(int i = 0; i != nonCriticalReactions.size(); ++i)
			std::cout << __func__ << ":" << __LINE__ << ": " << nonCriticalReactions[i] << " " << firings(i) << std::endl;
    }

    SKIP

    boost::numeric::ublas::vector<double> nonCriticalDeltas = boost::numeric::ublas::prod(stoichiometric, firings);
    deltas += nonCriticalDeltas;

    TIMER

    if(VERBOSE) {
        std::cout << __func__ << ":" << __LINE__ << ": deltas:";
        for(int i = 0; i != m.getNet().getNet().numPlaces(); ++i)
			std::cout << " " << deltas(i);
    	std::cout << std::endl;
    }

    SKIP

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

    TIMER

    PRINT_TIMINGS

    return {action, tau, true};
}

// ==============================================================================================

DrawMassActionEulerMaruyamaFunction::DrawMassActionEulerMaruyamaFunction(
	const lib::DG::Hyper &dg,
	std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> inputRate,
	std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> reactionRate,
	std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> outputRate,
	double tau)
	: dg(dg), inputRate(inputRate), reactionRate(reactionRate), outputRate(outputRate), tau(tau) {
	syncSize();
}

void DrawMassActionEulerMaruyamaFunction::syncSize() {
	const auto &g = dg.getGraph();
	const auto n = num_vertices(g);
	cachedInputRates.resize(n, -1.0);
	cachedRates.resize(n, -1.0);
}

std::tuple<Action, double, bool> DrawMassActionEulerMaruyamaFunction::draw(const Marking &m) {
	return draw_v0(m);
}

void DrawMassActionEulerMaruyamaFunction::syncState(const Marking &m) {
	const auto numPlaces = m.getNet().getNet().numPlaces();
	const auto oldSize = state.size();
	state.resize(numPlaces, true);
	if(!stateInitialised) {
		for(const auto v: asRange(vertices(dg.getGraph()))) {
			if(dg.getGraph()[v].kind != lib::DG::HyperVertexKind::Vertex) continue;
			const auto place = m.getNet().getPlace(v);
			state(place.getId()) = static_cast<double>(m.getMarking()[place]);
		}
		stateInitialised = true;
		return;
	}
	for(auto i = oldSize; i < numPlaces; ++i)
		state(i) = 0.0;
}

boost::numeric::ublas::vector<double> DrawMassActionEulerMaruyamaFunction::propensities(
		const Marking &m, const std::vector<int> &reactions) {
	const auto &dgGraph = dg.getGraph();
	boost::numeric::ublas::vector<double> results(reactions.size(), 0.0);
	for(std::size_t i = 0; i < reactions.size(); ++i) {
		const int idx = reactions[i];
		if(idx >= 0) {
			const auto v = vertices(dgGraph).first[idx];
			double reactionPropensity = 0.0;
			if(dgGraph[v].kind == lib::DG::HyperVertexKind::Edge) {
				reactionPropensity = realReactionPropensity(v, m, state);
				if(reactionPropensity == 0.0) continue;
			} else {
				const auto place = m.getNet().getPlace(v);
				if(state(place.getId()) <= 0.0) continue;
			}
			double r = cachedRates[idx];
			if(r < 0) {
				if(dgGraph[v].kind == lib::DG::HyperVertexKind::Edge && reactionRate) {
					bool cache;
					std::tie(r, cache) = reactionRate(dg, v);
					if(cache) cachedRates[idx] = r;
				} else if(dgGraph[v].kind == lib::DG::HyperVertexKind::Vertex && outputRate) {
					bool cache;
					std::tie(r, cache) = outputRate(dg, v);
					if(cache) cachedRates[idx] = r;
				} else {
					r = dgGraph[v].kind == lib::DG::HyperVertexKind::Edge ? 1.0 : 0.0;
					cachedRates[idx] = r;
				}
			}
			if(dgGraph[v].kind == lib::DG::HyperVertexKind::Edge) {
				results(i) = r * reactionPropensity;
			} else {
				const auto place = m.getNet().getPlace(v);
				results(i) = r * std::max(0.0, state(place.getId()));
			}
		} else {
			const auto rateIdx = -idx - 1;
			const auto v = vertices(dgGraph).first[rateIdx];
			double r = cachedInputRates[rateIdx];
			if(r < 0) {
				if(inputRate) {
					bool cache;
					std::tie(r, cache) = inputRate(dg, v);
					if(cache) cachedInputRates[rateIdx] = r;
				} else {
					cachedInputRates[rateIdx] = r = 0.0;
				}
			}
			results(i) = r;
		}
	}
	return results;
}

Action DrawMassActionEulerMaruyamaFunction::makeSyncAction(const Marking &m) const {
	std::vector<std::pair<lib::DG::HyperVertex, int>> updates;
	for(const auto v: asRange(vertices(dg.getGraph()))) {
		if(dg.getGraph()[v].kind != lib::DG::HyperVertexKind::Vertex) continue;
		const auto place = m.getNet().getPlace(v);
		const int desired = static_cast<int>(std::lround(std::max(0.0, state(place.getId()))));
		const int delta = desired - m.getMarking()[place];
		if(delta != 0) updates.emplace_back(v, delta);
	}
	return UpdateAction{std::move(updates)};
}

std::tuple<Action, double, bool> DrawMassActionEulerMaruyamaFunction::draw_v0(const Marking &m) {
	const auto &dgGraph = dg.getGraph();
	syncState(m);

	// idx of the reactions
	std::vector<int> reactions;
	for(const auto v: asRange(vertices(dgGraph))) {
		const auto idx = get(boost::vertex_index_t(), dgGraph, v);
		if(dgGraph[v].kind == lib::DG::HyperVertexKind::Edge) {
			reactions.emplace_back(idx);
		} else if(dgGraph[v].kind == lib::DG::HyperVertexKind::Vertex) {
			reactions.emplace_back(idx);
			reactions.emplace_back(-idx - 1);
		}
	}
	auto propensities = this->propensities(m, reactions);
	if(propensities.size() == 0 || !hasPositiveEntry(propensities))
		return {{}, 0.0, false};

	boost::numeric::ublas::mapped_matrix<double> stoichiometric(m.getNet().getNet().numPlaces(), reactions.size());
	for(std::size_t reaction = 0; reaction < reactions.size(); ++reaction) {
		const auto idx = reactions[reaction];
		for(const auto &[place, w] : consumed(dg, m, idx))
			stoichiometric(place.getId(), reaction) -= w;
		for(const auto &[place, w] : produced(dg, m, idx))
			stoichiometric(place.getId(), reaction) += w;
	}

	boost::numeric::ublas::vector<double> drift = boost::numeric::ublas::prod(stoichiometric, propensities);

	boost::numeric::ublas::vector<double> wienerIncrement(propensities.size());
	for(int i = 0; i < wienerIncrement.size(); i++) {
		std::normal_distribution<> dist(0, 1);
		auto &rng = mod::lib::getRng();
		const double rnd = dist(rng);
		wienerIncrement(i) = std::sqrt(propensities(i)) * rnd;
	}
	boost::numeric::ublas::vector<double> diffusion = boost::numeric::ublas::prod(stoichiometric, wienerIncrement);

	boost::numeric::ublas::vector<double> deltas = tau * drift + std::sqrt(tau) * diffusion;
	state += deltas;
	for(std::size_t i = 0; i < state.size(); ++i)
		state(i) = std::max(0.0, state(i));

	return {makeSyncAction(m), tau, true};
}

// ==============================================================================================

DrawMassActionSKRockFunction::DrawMassActionSKRockFunction(
	const lib::DG::Hyper &dg,
	std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> inputRate,
	std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> reactionRate,
	std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> outputRate,
	double tau,
	int stages)
	: dg(dg), inputRate(inputRate), reactionRate(reactionRate), outputRate(outputRate), tau(tau), stages(stages) {
	syncSize();
}

void DrawMassActionSKRockFunction::syncSize() {
	const auto &g = dg.getGraph();
	const auto n = num_vertices(g);
	cachedInputRates.resize(n, -1.0);
	cachedRates.resize(n, -1.0);
}

std::tuple<Action, double, bool> DrawMassActionSKRockFunction::draw(const Marking &m) {
	return draw_v0(m);
}

void DrawMassActionSKRockFunction::syncState(const Marking &m) {
	const auto numPlaces = m.getNet().getNet().numPlaces();
	const auto oldSize = state.size();
	state.resize(numPlaces, true);
	if(!stateInitialised) {
		for(const auto v: asRange(vertices(dg.getGraph()))) {
			if(dg.getGraph()[v].kind != lib::DG::HyperVertexKind::Vertex) continue;
			const auto place = m.getNet().getPlace(v);
			state(place.getId()) = static_cast<double>(m.getMarking()[place]);
		}
		stateInitialised = true;
		return;
	}
	for(auto i = oldSize; i < numPlaces; ++i)
		state(i) = 0.0;
}

Action DrawMassActionSKRockFunction::makeSyncAction(const Marking &m) const {
	std::vector<std::pair<lib::DG::HyperVertex, int>> updates;
	for(const auto v: asRange(vertices(dg.getGraph()))) {
		if(dg.getGraph()[v].kind != lib::DG::HyperVertexKind::Vertex) continue;
		const auto place = m.getNet().getPlace(v);
		const int desired = static_cast<int>(std::lround(std::max(0.0, state(place.getId()))));
		const int delta = desired - m.getMarking()[place];
		if(delta != 0) updates.emplace_back(v, delta);
	}
	return UpdateAction{std::move(updates)};
}

double DrawMassActionSKRockFunction::reactionPropensityWithDeltas(
    const Marking &m,
    lib::DG::HyperVertex e,
    boost::numeric::ublas::vector<double> deltas) {
	return realReactionPropensity(e, m, state, &deltas);
}

boost::numeric::ublas::vector<double> DrawMassActionSKRockFunction::propensitiesWithDeltas(
    const Marking &m,
    const std::vector<int> &reactions,
    boost::numeric::ublas::vector<double> deltas) {
    const auto &dgGraph = dg.getGraph();
    boost::numeric::ublas::vector<double> results(reactions.size(), 0.0);

    for(std::size_t i = 0; i < reactions.size(); i++) {
        const int idx = reactions[i];

        if(idx >= 0) {
            const auto v = vertices(dgGraph).first[idx];
            double reactionPropensity = 0.0;
            if(dgGraph[v].kind == lib::DG::HyperVertexKind::Edge) {
                reactionPropensity = reactionPropensityWithDeltas(m, v, deltas);
                if(reactionPropensity == 0.0) continue;
            } else {
                const auto place = m.getNet().getPlace(v);
                if(state(place.getId()) + deltas(place.getId()) <= 0.0) continue;
            }

            double r = cachedRates[idx];
            if(r < 0) {
                if(dgGraph[v].kind == lib::DG::HyperVertexKind::Edge && reactionRate) {
                    bool cache;
                    std::tie(r, cache) = reactionRate(dg, v);
                    if(cache) cachedRates[idx] = r;
                } else if(dgGraph[v].kind == lib::DG::HyperVertexKind::Vertex && outputRate) {
                    bool cache;
                    std::tie(r, cache) = outputRate(dg, v);
                    if(cache) cachedRates[idx] = r;
                } else {
                    r = dgGraph[v].kind == lib::DG::HyperVertexKind::Edge ? 1.0 : 0.0;
                    cachedRates[idx] = r;
                }
            }

            if(dgGraph[v].kind == lib::DG::HyperVertexKind::Edge)
                results(i) = r * reactionPropensity;
            else {
                const auto place = m.getNet().getPlace(v);
                results(i) = r * std::max(0.0, state(place.getId()) + deltas(place.getId()));
            }
        } else {
            const auto v = vertices(dgGraph).first[-idx - 1];

            double r = cachedInputRates[-idx - 1];
            if(r < 0) {
                if(inputRate) {
                    bool cache;
                    std::tie(r, cache) = inputRate(dg, v);
                    if(cache) cachedInputRates[-idx - 1] = r;
                } else {
                    cachedInputRates[-idx - 1] = r = 0.0;
                }
            }

            results(i) = r;
        }
	}

	return results;
}

// computes the function f when considering the deltas from the current marking
boost::numeric::ublas::vector<double> DrawMassActionSKRockFunction::f(
    const Marking &m,
    const boost::numeric::ublas::mapped_matrix<double> &stoichiometric,
    const std::vector<int> &reactions,
    boost::numeric::ublas::vector<double> deltas) {
    return boost::numeric::ublas::prod(stoichiometric, propensitiesWithDeltas(m, reactions, deltas));
}

std::tuple<Action, double, bool> DrawMassActionSKRockFunction::draw_v0(const Marking &m) {
    const auto &dgGraph = dg.getGraph();
	syncState(m);
	if(stages < 1)
		return {{}, 0.0, false};

	// idx of the reactions
	std::vector<int> reactions;
	for(const auto v: asRange(vertices(dgGraph))) {
		const auto idx = get(boost::vertex_index_t(), dgGraph, v);
		if(dgGraph[v].kind == lib::DG::HyperVertexKind::Edge) {
			reactions.emplace_back(idx);
		} else if(dgGraph[v].kind == lib::DG::HyperVertexKind::Vertex) {
			reactions.emplace_back(idx);
			reactions.emplace_back(-idx - 1);
		}
	}
	boost::numeric::ublas::vector<double> zeroDeltas(m.getNet().getNet().numPlaces(), 0.0);
	auto propensities = propensitiesWithDeltas(m, reactions, zeroDeltas);
	if(propensities.size() == 0 || !hasPositiveEntry(propensities))
		return {{}, 0.0, false};

	boost::numeric::ublas::mapped_matrix<double> stoichiometric(m.getNet().getNet().numPlaces(), reactions.size());
	for(std::size_t reaction = 0; reaction < reactions.size(); ++reaction) {
		const auto idx = reactions[reaction];
		for(const auto &[place, w] : consumed(dg, m, idx))
			stoichiometric(place.getId(), reaction) -= w;
		for(const auto &[place, w] : produced(dg, m, idx))
			stoichiometric(place.getId(), reaction) += w;
	}

    double eta = 0.05;
    double omega0 = 1 + eta / (static_cast<double>(stages) * static_cast<double>(stages));

    // precompute the values of the chebyshev polynomials evaluated at omega0
    std::vector<double> chebyshev(stages+1);
    chebyshev[0] = 1;
    chebyshev[1] = omega0;
    for(int i = 2; i <= stages; i++)
        chebyshev[i] = 2*omega0*chebyshev[i-1]-chebyshev[i-2];

    // precompute the values of the derivatives of the chebyshev polynomials evaluated at omega0
    std::vector<double> dchebyshev(stages+1);
    dchebyshev[0] = 0;
    dchebyshev[1] = 1;
    for(int i = 2; i <= stages; i++)
        dchebyshev[i] = 2*omega0*dchebyshev[i-1]+2*chebyshev[i-1]-dchebyshev[i-2];

	double omega1 = chebyshev[stages] / dchebyshev[stages];

    // precompute the values for mu, nu and kappa
    // TODO: optimize kappa away (can be computed from nu)
	std::vector<double> mus(stages), nus(stages), kappas(stages);
	mus[0] = omega1 / omega0;
	nus[0] = static_cast<double>(stages) * omega1 / 2;
    kappas[0] = static_cast<double>(stages) * omega1 / omega0;
	for(int i = 2; i <= stages; i++) {
		mus[i-1] = 2 * omega1 * chebyshev[i-1] / chebyshev[i];
		nus[i-1] = 2 * omega0 * chebyshev[i-1] / chebyshev[i];
		kappas[i-1] = - chebyshev[i-2] / chebyshev[i];
	}

	// compute the diffusion term Q
	boost::numeric::ublas::vector<double> wienerIncrement(propensities.size());
	for(int i = 0; i < wienerIncrement.size(); i++) {
		std::normal_distribution<> dist(0, 1);
		auto &rng = mod::lib::getRng();
		const double rnd = dist(rng);
		wienerIncrement(i) = std::sqrt(propensities(i)) * rnd;
	}
	boost::numeric::ublas::vector<double> diffusion = boost::numeric::ublas::prod(stoichiometric, wienerIncrement);
	boost::numeric::ublas::vector<double> Q = std::sqrt(tau) * diffusion;

    // compute the Ks
    std::vector<boost::numeric::ublas::vector<double>> Ks(stages+1);
    Ks[0] = boost::numeric::ublas::vector<double>(stoichiometric.size1(), 0.0);
    Ks[1] = mus[0] * tau * f(m, stoichiometric, reactions, nus[0] * Q) + kappas[0] * Q;
    for(int i = 2; i <= stages; i++) {
        boost::numeric::ublas::vector<double> fResult = f(m, stoichiometric, reactions, Ks[i-1]);
        Ks[i] = mus[i-1] * tau * fResult + nus[i-1] * Ks[i-1] + kappas[i-1] * Ks[i-2];
    }

    boost::numeric::ublas::vector<double> deltas = Ks[stages];
	state += deltas;
	for(std::size_t i = 0; i < state.size(); ++i)
		state(i) = std::max(0.0, state(i));

	return {makeSyncAction(m), tau, true};
}

// ==============================================================================================

void Simulator::doIteration() {
	++iteration;
}

} // namespace mod::lib::Causality
