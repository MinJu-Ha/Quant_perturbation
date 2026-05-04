from .base import PerturbationAdapter


class ScGPTAdapter(PerturbationAdapter):
    model_name = "scgpt"

    def run(self):
        try:
            import scgpt  # noqa: F401
        except ImportError:
            return self.result(
                "failed",
                "scgpt is not installed. Use mode: external_command or install/configure scGPT.",
            )

        required = ["checkpoint_dir", "vocab_path"]
        missing = [key for key in required if key not in self.runner_config]
        if missing:
            return self.result(
                "failed",
                f"scGPT python_api adapter needs these config keys: {', '.join(missing)}.",
            )

        return self.result(
            "failed",
            "scGPT package is importable, but the project-specific AnnData/checkpoint "
            "perturbation wrapper has not been configured yet.",
        )
