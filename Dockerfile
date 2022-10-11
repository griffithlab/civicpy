FROM python:3.7

RUN --mount=type=secret,id=civicpy_version \
  civicpy_version="$(cat /run/secrets/civicpy_version)" \
  && pip install civicpy==${civicpy_version}
