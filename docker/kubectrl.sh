kubectl create -f ./testbwa.yaml. #create the pod
kubectl get pods --show-all #look up which pod is the testbwa, get the name of pod

kubectl logs testbwa-#id #output the logs of the pod
kubectl delete jobs/testbwa #delete the pod

