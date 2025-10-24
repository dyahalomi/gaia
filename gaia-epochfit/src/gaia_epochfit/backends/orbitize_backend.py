class OrbitizeBackend:
    def __init__(self):
        try:
            import orbitize  # noqa: F401
        except Exception as e:
            raise ImportError("Install orbitize to use Orbitize backend") from e
    # Placeholder: Orbitize specializes in relative astrometry & DI.
