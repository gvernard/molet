// 'development' installs all 3rd party libraries.
// 'production' is the main image that installs the latest molet version from github.
// 'min' is the minified image of 'production' to save some space.
// If I modify vkl_lib or gerlumphpp then I need to rerun everything.
// If I modify only molet, then I just need to rerun 'production' and 'min'
//   Run the following in order:

docker build -t gvernard/molet:development . --network=host
docker build -t gvernard/molet:production . --network=host
docker build -t gvernard/molet:min . --network=host
docker push gvernard/molet:min
docker push gvernard/molet:production
docker push gvernard/molet:development
