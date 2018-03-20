# 峰科FPGA加速基因分析镜像设置说明

峰科基因镜像包括如下组件：
- 峰科基因分析软件*fcs-genome*，安装在`FALCON_DIR=/usr/local/falcon`
- 基于华为云的自动配置脚本，安装在`/root/setup-scripts/setup.sh`
- 峰科软件所需系统组件，相关环境变量配置存储在`/etc/profile.d/falcon.sh`

# 配置步骤
1. 申请一台32核ECS云服务器（例如: m2.8xlarge.8），和一台FP1云服务器。两台云服务器均选择 峰科基因镜像，且需要放在同一个可用区。
1. （可选）如需本地存储，请在ECS云服务器上挂在云盘，FP1云服务器不需挂载。
1. 配置之前，请记录如下信息：
    - FP1云服务器的私有IP（例如：192.168.x.x)
    - 云盘挂载位置（例如：/dev/vdb)
    - （可选）参考基因组位置
1. 登录ECS云服务器，执行自动配置脚本：
    ```
    > cd /root/setup-scripts
    > ./setup.sh
    ```    
    此配置脚本会执行如下几步：
    1. 配置SSH，需要输入FP1的私有IP。假设IP是*192.168.1.2*，配置界面如下：
        ```
        Please enter the private IP address of the FP1 node
        (e.g. 192.168.x.x): 192.168.1.2
        ```
    1. 配置FP1云服务器，无需用户输入
    1. 配置工作文件夹，此文件夹存储参考基因组和输入/输出结果，且需要对ECS云服务器 和FP1服务器共享。脚本提示如下：
        ```
        Setting up working dir...
        If you already have the working directory ready please enter
        the dir path, otherwise, please enter 'continue' or 'c':
        ```
    如果用户已提前配置好共享存储，此处直接输入共享存储地址。否则输入`continue`或者`c`。 下一步按照提示输入本地存储挂载位置和工作目录地址，以下假设挂载位置为`/dev/vdb`， 工作文件夹为`WORK_DIR=/local`：
        ```
        If you already have the working directory ready please enter
        the dir path, otherwise, please enter 'continue' or 'c': c

        Please enter the storage device: /dev/vdb
        Please enter the path to the work dir: /local
        ```
    1. （可选）预处理参考基因组，自动化脚本可以初始化参考基因组，如需跳过此步可直接回车不输入参考基因位置：
        ```
        Preparing reference genome...
        Please enter the path to the reference genome, or leave it blank to skip this step: /local/ref/human_g1k_v37.fasta
        ```
    1. 配置完毕会出现如下信息：
        ```
        Configuration is successful.
        ```
    配置成功后，查看`/usr/local/falcon/fcs-genome.conf`配置文件将会有如下两行：
        ```
        hosts = localhost,192.168.1.176
        temp_dir = /local/temp
        ```

# 运行前准备
1. 确保工作目录有足够大的存储空间，并且确包参考基因组和输入文件存储在此工作目录中。一般情况下工作目录下的可用空间应为输入*FASTQ.GZ*文件的3到5倍。
1. 确保峰科软件License成功配置在环境变量`$LM_LICENSE_PATH`中。此环境变量既可以指向文件，也可以是远程地址：`<port>@<hostname>`。如果License配置错误或者License已过期，将会出现 如下提示：
    ```
    [fcs-genome] ERROR: Cannot connect to the license server: -15
    [fcs-genome] ERROR: Please contact support@falcon-computing.com for details.
    ```
    如需获取最新License请联系峰科销售，或者发邮件至support@falcon-computing.com
1. （可选）预处理参考基因一遍最大化峰科基因加速软件的性能。普通参考基因组未经预处理也可以使用，所以此步可以跳过。处理过的参考基因组也与标准文件匹配，不影响其他基因分析软件如*BWA*, *GATK*, *Picard*等的使用。预处理可以通过如下命令执行：
    ```
    > $FALCON_DIR/prepare-ref.sh <path-to-fasta>
    ```
1. （可选）GATK数据库如*dbSNP*等如果没有索引文件，会影响后续基因分析流程的性能。所以建议提前生成索引文件，可以使用标准*GATK*或者[*samtools distribution*](http://www.htslib.org/download/)中的*tabix*工具。

# 运行软件
本镜像提供人重基因组分析示例脚本，位置在`$FALCON_DIR/example-wgs-germline.sh`。 具体使用方法请参考峰科基因软件`Falcon Genome`的使用文档（英文）。
