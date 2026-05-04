import subprocess

from .base import PerturbationAdapter


class ExternalCommandAdapter(PerturbationAdapter):
    model_name = "external"

    def run(self):
        command = self.runner_config.get("command")
        if not command:
            return self.result("failed", "No external command configured.")

        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0:
            return self.result("failed", f"Command exited with code {completed.returncode}.")

        if not self.output_path.exists():
            return self.result("failed", f"Command finished but did not write {self.output_path}.")

        return self.result("completed", f"Wrote {self.output_path}.")
