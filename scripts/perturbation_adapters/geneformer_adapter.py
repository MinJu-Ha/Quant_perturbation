from .base import PerturbationAdapter


class GeneformerAdapter(PerturbationAdapter):
    model_name = "geneformer"

    def run(self):
        try:
            from geneformer import InSilicoPerturber  # noqa: F401
        except ImportError:
            return self.result(
                "failed",
                "geneformer is not installed. Use mode: external_command or install/configure Geneformer.",
            )

        required = ["model_directory", "tokenized_dataset_path"]
        missing = [key for key in required if key not in self.runner_config]
        if missing:
            return self.result(
                "failed",
                f"Geneformer python_api adapter needs these config keys: {', '.join(missing)}.",
            )

        return self.result(
            "failed",
            "Geneformer predicts embedding/state shifts rather than a native genes x cells "
            "expression matrix. Add a projection/export step before scoring with this evaluator.",
        )
