"""
Base classes for pipeline implementation
"""

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Any, Optional, Union, List
from enum import Enum

class StepStatus(Enum):
    """Status of a pipeline step"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"

class PipelineStep(ABC):
    """Base class for a single pipeline step"""
    
    def __init__(self, name: str, workdir: Path, config: Dict[str, Any] = None):
        """
        Initialize pipeline step
        
        Args:
            name: Step name
            workdir: Working directory for this step
            config: Configuration dictionary
        """
        self.name = name
        self.workdir = Path(workdir) / name
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.config = config or {}
        self.logger = logging.getLogger(f"{__name__}.{name}")
        self.status = StepStatus.PENDING
        self.outputs: Dict[str, Path] = {}
        
    @abstractmethod
    def _execute(self, **kwargs) -> Dict[str, Any]:
        """
        Execute step logic - to be implemented by subclasses
        
        Returns:
            Dictionary with step results
        """
        pass
    
    def run(self, **kwargs) -> Dict[str, Any]:
        """
        Run the step
        
        Returns:
            Dictionary with step results including status and outputs
        """
        self.status = StepStatus.RUNNING
        self.logger.info(f"Starting step: {self.name}")
        
        try:
            result = self._execute(**kwargs)
            self.status = StepStatus.COMPLETED
            self.logger.info(f"Step {self.name} completed")
            
            # Add metadata to result
            result['status'] = self.status.value
            result['step_name'] = self.name
            result['workdir'] = str(self.workdir)
            result['outputs'] = {k: str(v) for k, v in self.outputs.items()}
            
            return result
            
        except Exception as e:
            self.status = StepStatus.FAILED
            self.logger.error(f"Step {self.name} failed: {e}")
            raise
    
    def get_output(self, key: str) -> Optional[Path]:
        """
        Get output file by key
        
        Args:
            key: Output identifier
            
        Returns:
            Path to output file or None
        """
        return self.outputs.get(key)


class Pipeline(ABC):
    """Base class for a complete pipeline with multiple steps"""
    
    def __init__(self, workdir: Union[str, Path], config: Dict[str, Any] = None):
        """
        Initialize pipeline
        
        Args:
            workdir: Working directory for the pipeline
            config: Configuration dictionary
        """
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.config = config or {}
        self.logger = logging.getLogger(__name__)
        self.steps: Dict[str, PipelineStep] = {}
        self.results: Dict[str, Any] = {}
        
    def add_step(self, name: str, step: PipelineStep):
        """
        Add a step to the pipeline
        
        Args:
            name: Step name
            step: PipelineStep instance
        """
        self.steps[name] = step
        
    def run_step(self, name: str, **kwargs) -> Dict[str, Any]:
        """
        Run a single step by name
        
        Args:
            name: Step name
            **kwargs: Arguments for the step
            
        Returns:
            Step results
        """
        if name not in self.steps:
            raise KeyError(f"Step not found: {name}")
        
        self.logger.info(f"Running step: {name}")
        result = self.steps[name].run(**kwargs)
        self.results[name] = result
        return result
    
    def run(self, steps: List[str] = None, **kwargs) -> Dict[str, Any]:
        """
        Run specified steps or all steps
        
        Args:
            steps: List of step names to run (None for all)
            **kwargs: Arguments for steps
            
        Returns:
            Dictionary with all step results
        """
        steps_to_run = steps or list(self.steps.keys())
        for name in steps_to_run:
            self.run_step(name, **kwargs)
        return self.results


