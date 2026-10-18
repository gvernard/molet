// 'development' installs all 3rd party libraries.
// 'production' is the main image that installs the latest molet version from github.
// If I modify vkl_lib or gerlumphpp then I need to rerun everything.
// If I modify only molet, then I just need to rerun 'development' and 'production'.
//   Run the following in order:

docker build -t gvernard/molet:development . --network=host
docker build -t gvernard/molet:production . --network=host
docker push gvernard/molet:production
docker push gvernard/molet:development
