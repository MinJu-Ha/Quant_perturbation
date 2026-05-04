from dataclasses import dataclass
from pathlib import Path


@dataclass
class AdapterResult:
    model: str
    output_path: Path
    output_format: str
    status: str
    message: str


class PerturbationAdapter:
    model_name = "base"

    def __init__(self, config, runner_config):
        self.config = config
        self.runner_config = runner_config

    @property
    def output_path(self):
        return Path(self.runner_config["output_path"])

    @property
    def output_format(self):
        return self.runner_config["output_format"]

    def run(self):
        raise NotImplementedError

    def result(self, status, message):
        return AdapterResult(
            model=self.model_name,
            output_path=self.output_path,
            output_format=self.output_format,
            status=status,
            message=message,
        )
