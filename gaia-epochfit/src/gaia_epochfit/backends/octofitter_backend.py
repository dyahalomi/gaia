class OctofitterBackend:
    def __init__(self):
        try:
            import octofitterpy  # noqa: F401
        except Exception as e:
            raise ImportError("Install octofitterpy to use Octofitter backend") from e
    # Sketch only; use octofitterpy to pass absolute astrometry or acceleration terms.
