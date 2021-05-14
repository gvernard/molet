#!/bin/bash
cp Dockerfile_full Dockerfile
docker build -t gvernard/molet:development . --network=host
cp Dockerfile_production Dockerfile
docker build -t gvernard/molet:production . --network=host
cp Dockerfile_min Dockerfile
docker build -t gvernard/molet:min . --network=host
docker push gvernard/molet:min
docker push gvernard/molet:production
docker push gvernard/molet:development
rm Dockerfile
