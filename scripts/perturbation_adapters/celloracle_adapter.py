from .base import PerturbationAdapter


class CellOracleAdapter(PerturbationAdapter):
    model_name = "celloracle"

    def run(self):
        try:
            import celloracle  # noqa: F401
        except ImportError:
            return self.result(
                "failed",
                "celloracle is not installed. Use mode: external_command or install/configure CellOracle.",
            )

        return self.result(
            "failed",
            "CellOracle python_api adapter needs a fitted oracle/GRN object path. "
            "Use an external runner for now, or add oracle_path to model_runners.celloracle.",
        )
