#include "Stochsim.hpp"

#include <mod/Error.hpp>
#include <mod/causality/Petri.hpp>
#include <mod/causality/EventTrace.hpp>
#include <mod/lib/Causality/Stochsim.hpp>

namespace mod::causality {

struct DrawMassActionFunction::Pimpl {
	std::shared_ptr<dg::DG> dg_;
	lib::Causality::DrawMassActionFunction m;
};

DrawMassActionFunction::DrawMassActionFunction(std::shared_ptr<dg::DG> dg_,
                                               std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate,
                                               std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate,
                                               std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate) {
	if(!dg_) throw LogicError("The derivation graph is a null pointer.");
	if(!dg_->hasActiveBuilder() && !dg_->isLocked())
		throw LogicError("The DG neither has an active builder nor is locked yet.");

	using F = std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)>;
	F inputRateInner, reactionRateInner, outputRateInner;
	if(inputRate) {
		inputRateInner = [inputRate](const lib::DG::Hyper &dgHyper,
		                             const lib::DG::HyperVertex v) -> std::pair<double, bool> {
			return inputRate(dgHyper.getInterfaceVertex(v));
		};
	}
	if(reactionRate) {
		reactionRateInner = [reactionRate](const lib::DG::Hyper &dgHyper,
		                                   const lib::DG::HyperVertex e) -> std::pair<double, bool> {
			return reactionRate(dgHyper.getInterfaceEdge(e));
		};
	}
	if(outputRate) {
		outputRateInner = [outputRate](const lib::DG::Hyper &dgHyper,
		                               const lib::DG::HyperVertex v) -> std::pair<double, bool> {
			return outputRate(dgHyper.getInterfaceVertex(v));
		};
	}
	p.reset(new Pimpl{
		dg_, lib::Causality::DrawMassActionFunction(
				dg_->getHyper(), inputRateInner, reactionRateInner, outputRateInner)
	});
}

DrawMassActionFunction::~DrawMassActionFunction() = default;
DrawMassActionFunction::DrawMassActionFunction(DrawMassActionFunction &&) = default;
DrawMassActionFunction &DrawMassActionFunction::operator=(DrawMassActionFunction &&) = default;

DrawMassActionFunction::DrawMassActionFunction(const DrawMassActionFunction &other) {
	p.reset(new Pimpl(*other.p));
}

DrawMassActionFunction &DrawMassActionFunction::operator=(const DrawMassActionFunction &other) {
	if(&other != this)
		p.reset(new Pimpl(*other.p));
	return *this;
}

void DrawMassActionFunction::syncSize() {
	p->m.syncSize();
}

std::tuple<std::optional<Action>, double, bool> DrawMassActionFunction::draw(const Marking &m) {
	if(m.getNet()->getDG() != p->dg_) throw LogicError("The marking is not on the underlying derivation graph.");
	const auto [actionInner, total, isTimeIncrement] = p->m.draw(m.getMarking());
	if(total == 0) return {std::nullopt, 0, false};
	struct Convert {
		Action operator()(lib::Causality::EdgeAction a) const {
			return EdgeAction(dgHyper.getInterfaceEdge(a.e));
		}

		Action operator()(lib::Causality::InputAction a) const {
			return InputAction(dgHyper.getInterfaceVertex(a.v));
		}

		Action operator()(lib::Causality::OutputAction a) const {
			return OutputAction(dgHyper.getInterfaceVertex(a.v));
		}

		Action operator()(lib::Causality::UpdateAction a) const {
            std::vector<std::pair<dg::DG::Vertex, int>> updates;
            updates.reserve(a.updates.size());
            for (const auto &[v, c] : a.updates)
                updates.emplace_back(dgHyper.getInterfaceVertex(v), c);
            return UpdateAction(std::move(updates));
        }
	public:
		const lib::DG::Hyper &dgHyper;
	};
	return {std::visit(Convert{p->dg_->getHyper()}, actionInner), total, isTimeIncrement};
}

// =============================================================================================================

struct DrawMassActionTauLeapingFunction::Pimpl {
	std::shared_ptr<dg::DG> dg_;
	lib::Causality::DrawMassActionTauLeapingFunction m;
};

DrawMassActionTauLeapingFunction::DrawMassActionTauLeapingFunction(std::shared_ptr<dg::DG> dg_,
        std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate,
        std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate,
        std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate,
        int dc, double epsilon) {
	if(!dg_) throw LogicError("The derivation graph is a null pointer.");
	if(!dg_->hasActiveBuilder() && !dg_->isLocked())
		throw LogicError("The DG neither has an active builder nor is locked yet.");

	using F = std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)>;
	F inputRateInner, reactionRateInner, outputRateInner;
	if(inputRate) {
		inputRateInner = [inputRate](const lib::DG::Hyper &dgHyper,
		                             const lib::DG::HyperVertex v) -> std::pair<double, bool> {
			return inputRate(dgHyper.getInterfaceVertex(v));
		};
	}
	if(reactionRate) {
		reactionRateInner = [reactionRate](const lib::DG::Hyper &dgHyper,
		                                   const lib::DG::HyperVertex e) -> std::pair<double, bool> {
			return reactionRate(dgHyper.getInterfaceEdge(e));
		};
	}
	if(outputRate) {
		outputRateInner = [outputRate](const lib::DG::Hyper &dgHyper,
		                               const lib::DG::HyperVertex v) -> std::pair<double, bool> {
			return outputRate(dgHyper.getInterfaceVertex(v));
		};
	}
	p.reset(new Pimpl{
		dg_, lib::Causality::DrawMassActionTauLeapingFunction(
				dg_->getHyper(), inputRateInner, reactionRateInner, outputRateInner, dc, epsilon)
	});
}

DrawMassActionTauLeapingFunction::~DrawMassActionTauLeapingFunction() = default;
DrawMassActionTauLeapingFunction::DrawMassActionTauLeapingFunction(DrawMassActionTauLeapingFunction &&) = default;
DrawMassActionTauLeapingFunction &DrawMassActionTauLeapingFunction::operator=(DrawMassActionTauLeapingFunction &&) = default;

DrawMassActionTauLeapingFunction::DrawMassActionTauLeapingFunction(const DrawMassActionTauLeapingFunction &other) {
	p.reset(new Pimpl(*other.p));
}

DrawMassActionTauLeapingFunction &DrawMassActionTauLeapingFunction::operator=(const DrawMassActionTauLeapingFunction &other) {
	if(&other != this)
		p.reset(new Pimpl(*other.p));
	return *this;
}

void DrawMassActionTauLeapingFunction::syncSize() {
	p->m.syncSize();
}

std::tuple<std::optional<Action>, double, bool> DrawMassActionTauLeapingFunction::draw(const Marking &m) {
	if(m.getNet()->getDG() != p->dg_) throw LogicError("The marking is not on the underlying derivation graph.");
	const auto [actionInner, total, isTimeIncrement] = p->m.draw(m.getMarking());
	if(total == 0) return {std::nullopt, 0, false};
	struct Convert {
		Action operator()(lib::Causality::EdgeAction a) const {
			return EdgeAction(dgHyper.getInterfaceEdge(a.e));
		}

		Action operator()(lib::Causality::InputAction a) const {
			return InputAction(dgHyper.getInterfaceVertex(a.v));
		}

		Action operator()(lib::Causality::OutputAction a) const {
			return OutputAction(dgHyper.getInterfaceVertex(a.v));
		}

		Action operator()(lib::Causality::UpdateAction a) const {
            std::vector<std::pair<dg::DG::Vertex, int>> updates;
            updates.reserve(a.updates.size());
            for (const auto &[v, c] : a.updates)
                updates.emplace_back(dgHyper.getInterfaceVertex(v), c);
            return UpdateAction(std::move(updates));
        }
	public:
		const lib::DG::Hyper &dgHyper;
	};
	return {std::visit(Convert{p->dg_->getHyper()}, actionInner), total, isTimeIncrement};
}

// =============================================================================================================

struct DrawMassActionEulerMaruyamaFunction::Pimpl {
	std::shared_ptr<dg::DG> dg_;
	lib::Causality::DrawMassActionEulerMaruyamaFunction m;
};

DrawMassActionEulerMaruyamaFunction::DrawMassActionEulerMaruyamaFunction(std::shared_ptr<dg::DG> dg_,
        std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate,
        std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate,
        std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate,
        double tau) {
	if(!dg_) throw LogicError("The derivation graph is a null pointer.");
	if(!dg_->hasActiveBuilder() && !dg_->isLocked())
		throw LogicError("The DG neither has an active builder nor is locked yet.");

	using F = std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)>;
	F inputRateInner, reactionRateInner, outputRateInner;
	if(inputRate) {
		inputRateInner = [inputRate](const lib::DG::Hyper &dgHyper,
		                             const lib::DG::HyperVertex v) -> std::pair<double, bool> {
			return inputRate(dgHyper.getInterfaceVertex(v));
		};
	}
	if(reactionRate) {
		reactionRateInner = [reactionRate](const lib::DG::Hyper &dgHyper,
		                                   const lib::DG::HyperVertex e) -> std::pair<double, bool> {
			return reactionRate(dgHyper.getInterfaceEdge(e));
		};
	}
	if(outputRate) {
		outputRateInner = [outputRate](const lib::DG::Hyper &dgHyper,
		                               const lib::DG::HyperVertex v) -> std::pair<double, bool> {
			return outputRate(dgHyper.getInterfaceVertex(v));
		};
	}
	p.reset(new Pimpl{
		dg_, lib::Causality::DrawMassActionEulerMaruyamaFunction(
				dg_->getHyper(), inputRateInner, reactionRateInner, outputRateInner, tau)
	});
}

DrawMassActionEulerMaruyamaFunction::~DrawMassActionEulerMaruyamaFunction() = default;
DrawMassActionEulerMaruyamaFunction::DrawMassActionEulerMaruyamaFunction(DrawMassActionEulerMaruyamaFunction &&) = default;
DrawMassActionEulerMaruyamaFunction &DrawMassActionEulerMaruyamaFunction::operator=(DrawMassActionEulerMaruyamaFunction &&) = default;

DrawMassActionEulerMaruyamaFunction::DrawMassActionEulerMaruyamaFunction(const DrawMassActionEulerMaruyamaFunction &other) {
	p.reset(new Pimpl(*other.p));
}

DrawMassActionEulerMaruyamaFunction &DrawMassActionEulerMaruyamaFunction::operator=(const DrawMassActionEulerMaruyamaFunction &other) {
	if(&other != this)
		p.reset(new Pimpl(*other.p));
	return *this;
}

void DrawMassActionEulerMaruyamaFunction::syncSize() {
	p->m.syncSize();
}

std::tuple<std::optional<Action>, double, bool> DrawMassActionEulerMaruyamaFunction::draw(const Marking &m) {
	if(m.getNet()->getDG() != p->dg_) throw LogicError("The marking is not on the underlying derivation graph.");
	const auto [actionInner, total, isTimeIncrement] = p->m.draw(m.getMarking());
	if(total == 0) return {std::nullopt, 0, false};
	struct Convert {
		Action operator()(lib::Causality::EdgeAction a) const {
			return EdgeAction(dgHyper.getInterfaceEdge(a.e));
		}

		Action operator()(lib::Causality::InputAction a) const {
			return InputAction(dgHyper.getInterfaceVertex(a.v));
		}

		Action operator()(lib::Causality::OutputAction a) const {
			return OutputAction(dgHyper.getInterfaceVertex(a.v));
		}

		Action operator()(lib::Causality::UpdateAction a) const {
            std::vector<std::pair<dg::DG::Vertex, int>> updates;
            updates.reserve(a.updates.size());
            for (const auto &[v, c] : a.updates)
                updates.emplace_back(dgHyper.getInterfaceVertex(v), c);
            return UpdateAction(std::move(updates));
        }
	public:
		const lib::DG::Hyper &dgHyper;
	};
	return {std::visit(Convert{p->dg_->getHyper()}, actionInner), total, isTimeIncrement};
}

// =============================================================================================================

struct DrawMassActionSKRockFunction::Pimpl {
	std::shared_ptr<dg::DG> dg_;
	lib::Causality::DrawMassActionSKRockFunction m;
};

DrawMassActionSKRockFunction::DrawMassActionSKRockFunction(std::shared_ptr<dg::DG> dg_,
        std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate,
        std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate,
        std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate,
        double tau,
        int stages) {
	if(!dg_) throw LogicError("The derivation graph is a null pointer.");
	if(!dg_->hasActiveBuilder() && !dg_->isLocked())
		throw LogicError("The DG neither has an active builder nor is locked yet.");

	using F = std::function<std::pair<double, bool>(const lib::DG::Hyper &, lib::DG::HyperVertex)>;
	F inputRateInner, reactionRateInner, outputRateInner;
	if(inputRate) {
		inputRateInner = [inputRate](const lib::DG::Hyper &dgHyper,
		                             const lib::DG::HyperVertex v) -> std::pair<double, bool> {
			return inputRate(dgHyper.getInterfaceVertex(v));
		};
	}
	if(reactionRate) {
		reactionRateInner = [reactionRate](const lib::DG::Hyper &dgHyper,
		                                   const lib::DG::HyperVertex e) -> std::pair<double, bool> {
			return reactionRate(dgHyper.getInterfaceEdge(e));
		};
	}
	if(outputRate) {
		outputRateInner = [outputRate](const lib::DG::Hyper &dgHyper,
		                               const lib::DG::HyperVertex v) -> std::pair<double, bool> {
			return outputRate(dgHyper.getInterfaceVertex(v));
		};
	}
	p.reset(new Pimpl{
		dg_, lib::Causality::DrawMassActionSKRockFunction(
				dg_->getHyper(), inputRateInner, reactionRateInner, outputRateInner, tau, stages)
	});
}

DrawMassActionSKRockFunction::~DrawMassActionSKRockFunction() = default;
DrawMassActionSKRockFunction::DrawMassActionSKRockFunction(DrawMassActionSKRockFunction &&) = default;
DrawMassActionSKRockFunction &DrawMassActionSKRockFunction::operator=(DrawMassActionSKRockFunction &&) = default;

DrawMassActionSKRockFunction::DrawMassActionSKRockFunction(const DrawMassActionSKRockFunction &other) {
	p.reset(new Pimpl(*other.p));
}

DrawMassActionSKRockFunction &DrawMassActionSKRockFunction::operator=(const DrawMassActionSKRockFunction &other) {
	if(&other != this)
		p.reset(new Pimpl(*other.p));
	return *this;
}

void DrawMassActionSKRockFunction::syncSize() {
	p->m.syncSize();
}

std::tuple<std::optional<Action>, double, bool> DrawMassActionSKRockFunction::draw(const Marking &m) {
	if(m.getNet()->getDG() != p->dg_) throw LogicError("The marking is not on the underlying derivation graph.");
	const auto [actionInner, total, isTimeIncrement] = p->m.draw(m.getMarking());
	if(total == 0) return {std::nullopt, 0, false};
	struct Convert {
		Action operator()(lib::Causality::EdgeAction a) const {
			return EdgeAction(dgHyper.getInterfaceEdge(a.e));
		}

		Action operator()(lib::Causality::InputAction a) const {
			return InputAction(dgHyper.getInterfaceVertex(a.v));
		}

		Action operator()(lib::Causality::OutputAction a) const {
			return OutputAction(dgHyper.getInterfaceVertex(a.v));
		}

		Action operator()(lib::Causality::UpdateAction a) const {
            std::vector<std::pair<dg::DG::Vertex, int>> updates;
            updates.reserve(a.updates.size());
            for (const auto &[v, c] : a.updates)
                updates.emplace_back(dgHyper.getInterfaceVertex(v), c);
            return UpdateAction(std::move(updates));
        }
	public:
		const lib::DG::Hyper &dgHyper;
	};
	return {std::visit(Convert{p->dg_->getHyper()}, actionInner), total, isTimeIncrement};
}

// =============================================================================================================

struct SimulatorImpl::Pimpl {
	lib::Causality::Simulator sim;
};

SimulatorImpl::SimulatorImpl() : p(new Pimpl()) {}

SimulatorImpl::~SimulatorImpl() = default;

int SimulatorImpl::getIteration() const { return p->sim.getIteration(); }
double SimulatorImpl::getTime() const { return p->sim.getTime(); }
void SimulatorImpl::setTime_delete(double value) { p->sim.setTime_delete(value); }

void SimulatorImpl::doIteration() {
	p->sim.doIteration();
}

} // namespace mod::causality
