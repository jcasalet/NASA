#!/usr/bin/env bash
if docker build -t ah-causal-ensemble-docker .  ; then
  docker tag ah-causal-ensemble-docker gcr.io/fdl-astronaut-health/ah-causal-ensemble-docker
  docker push gcr.io/fdl-astronaut-health/ah-causal-ensemble-docker

  docker tag ah-causal-ensemble-docker registry.gitlab.com/frontierdevelopmentlab/astronaut-health/ah-causal-ensemble
  docker push registry.gitlab.com/frontierdevelopmentlab/astronaut-health/ah-causal-ensemble
fi
