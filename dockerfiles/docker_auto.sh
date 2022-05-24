#!/bin/bash
printf '%*s\n' "${COLUMNS:=$(tput cols)}" '' | tr ' ' -
echo "Docker tag: FULL"
docker rmi $(docker images 'gvernard/molet:development' -a -q)
cp Dockerfile_full Dockerfile
docker build --build-arg CACHE_DATE="$(date)" -t gvernard/molet:development . --network=host

printf '%*s\n' "${COLUMNS:=$(tput cols)}" '' | tr ' ' -
echo "Docker tag: PRODUCTION"
docker rmi $(docker images 'gvernard/molet:production' -a -q)
cp Dockerfile_prod Dockerfile
docker build --build-arg CACHE_DATE="$(date)" -t gvernard/molet:production . --network=host

printf '%*s\n' "${COLUMNS:=$(tput cols)}" '' | tr ' ' -
echo "Docker tag: MIN"
docker rmi $(docker images 'gvernard/molet:min' -a -q)
cp Dockerfile_min Dockerfile
docker build -t gvernard/molet:min . --network=host


printf '%*s\n' "${COLUMNS:=$(tput cols)}" '' | tr ' ' -
echo "Pushing docker images..."
echo "Pushing: min..."
docker push gvernard/molet:min
echo "DONE"

echo "Pushing: production..."
docker push gvernard/molet:production
echo "DONE"

echo "Pushing: development..."
docker push gvernard/molet:development
echo "DONE"

rm Dockerfile
echo "ALL DONE"
