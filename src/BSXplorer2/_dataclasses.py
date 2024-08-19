from dataclasses import dataclass, asdict

@dataclass
class DClass:
    def as_dict(self):
        return asdict(self)

@dataclass
class ChromosomeDClass(DClass):
    name: str
    start: int
    end: int
    length: int
