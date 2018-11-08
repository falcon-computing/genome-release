kubectl create -f ./fcs.yaml. #create the pod
kubectl get pods --show-all #look up which pod is the fcs, get the name of pod

kubectl logs fcs-#id #output the logs of the pod
kubectl delete jobs/fcs #delete the pod

