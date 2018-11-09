#build docker image from current folder, using ./Dockerfile
dockerName=$1
docker build -t $dockerName .

