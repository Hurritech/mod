#ifndef MOD_LIB_CAUSALITY_STOCHSIM_HPP
#define MOD_LIB_CAUSALITY_STOCHSIM_HPP

#include <mod/lib/Causality/EventTrace.hpp>

#include <functional>

namespace mod::lib::DG {
struct Hyper;
} // namespace mod::lib::DG
namespace mod::lib::Causality {

struct DrawMassActionFunction {
	DrawMassActionFunction(const lib::DG::Hyper &dg,
	                       std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)>
	                       inputRate,
	                       std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)>
	                       reactionRate,
	                       std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)>
	                       outputRate);
	void syncSize();
	// .second is 0.0 when no actions are possible
	std::pair<Action, double> draw(const Marking &m);
private:
	std::pair<Action, double> draw_v0(const Marking &m);
private:
	const lib::DG::Hyper &dg;
	const std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)>
			inputRate, reactionRate, outputRate;
	std::vector<double> cachedInputRates, cachedRates /* reaction and output */;
};

struct DrawMassActionTauLeapingFunction {
    DrawMassActionTauLeapingFunction(const lib::DG::Hyper &dg,
        std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> inputRate,
        std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> reactionRate,
        std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)> outputRate,
        int dc, double epsilon);
    void syncSize();
    std::pair<Action, double> draw(const Marking &m);
private:
    std::pair<Action, double> draw_v0(const Marking &m);
private:
    const lib::DG::Hyper &dg;
    const std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)>
            inputRate, reactionRate, outputRate;
    const int dc;
    const double epsilon;
    std::vector<double> cachedInputRates, cachedRates;
};

struct Simulator {
public:
	int getIteration() const { return iteration; }
	double getTime() const { return time; }
	void setTime_delete(double value) { time = value; }
public:
	void doIteration();
private:
	int iteration = 0;
	double time = 0;
};

} // namespace mod::lib::Causality

#endif // MOD_LIB_CAUSALITY_STOCHSIM_HPP
