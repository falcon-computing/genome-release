#! /usr/bin/python
import atexit
import os.path
import socket
import subprocess
import struct
import sys
import task_pb2
import time
from threading import *

def sock_send(s, msg):
  msg_len = struct.pack('<I', len(msg))
  msg_str = msg_len + msg
  s.send(msg_str.encode())

def sock_recv(s):
  l_str = s.recv(4).decode()
  l = struct.unpack('<I', l_str)
  ret = s.recv(l[0])
  return ret

if len(sys.argv) != 2: 
  print "Usage:", sys.argv[0], "acc_id"
  sys.exit(-1)

acc_id = sys.argv[1]
print "info: checking accelerator [" + acc_id + "]"

# check nam binary and conf
falcon_home = '/usr/local/falcon'
if 'FALCON_HOME' in os.environ:
  falcon_home = os.getenv('FALCON_HOME')
else:
  print "warning: missing env $FALCON_HOME, use default: /usr/local/falcon"

nam_path  = falcon_home + '/blaze/bin/nam'
conf_path = falcon_home + '/blaze/conf'

if not os.path.isfile(nam_path) or not os.path.isfile(conf_path): 
  print "error: missing blaze nam or configuration file"
  sys.exit(1)

# start blaze nam in the background
blaze_proc = subprocess.Popen([nam_path, conf_path])
print "info:", nam_path, conf_path


# kill nam at exit
def cleanup():
  p = blaze_proc
  timeout_sec = 5
  p_sec = 0
  for second in range(timeout_sec):
    if p.poll() == None:
      time.sleep(1)
      p_sec += 1

  if p_sec >= timeout_sec:
    p.kill() # supported from python 2.6

atexit.register(cleanup)

# generate a task request message
task = task_pb2.TaskMsg()
task.type = task_pb2.ACCREQUEST
task.acc_id = acc_id
task.app_id = 'test-blaze'
task_str = task.SerializeToString()

ret_str = ''

count = 0
while count < 5:
  try:
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect(('127.0.0.1', 1027))
  
    sock_send(s, task_str)
  
    ret_str = sock_recv(s)
  
    s.close()

    # socket connection was successful
    break
  
  except:
    # give blaze nam some time to start
    time.sleep(5)
    count += 1

if count == 4:
  print "error: failed to communicate with nam"
  sys.exit(-1)

# decode task message
ret = task_pb2.TaskMsg()
ret.ParseFromString(ret_str)
if ret.type != task_pb2.ACCGRANT:
  print "info: accelerator does not exist"
  sys.exit(1)
else:
  print "info: accelerator exists"
