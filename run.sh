docker run \
	--privileged \
	-v /data/isilon/users/berutti/Projects/Pipeline/Docker/Pipeline/configuration_test:/configuration \
	-v /data/isilon/users/berutti/Projects/Pipeline/Docker/Pipeline/testdata:/testdata \
	-v /data/mirror/goldenpath:/data/mirror/goldenpath \
	-h $(hostname)docker   \
	-it \
	test:0.1


# Docka socka
#	-v /var/run/docker.sock:/var/run/docker.sock \
