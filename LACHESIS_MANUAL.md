## LACHESIS 的安装流程
### 安装 samtools-0.1.18
去 SOURCEFOREGE 的 页面找到 [samtools-0.1.18](https://sourceforge.net/projects/samtools/files/samtools/0.1.18/) 的压缩包，下载下来。
```
tar xf samtools-0.1.18.tar.bz2
cd samtools-0.1.18
make -j 20
```

### 安装 boost 1.59.0
按照 [官方文档](https://www.boost.org/doc/libs/1_59_0/more/getting_started/unix-variants.html) 的指引进行安装：
```
wget http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.bz2
tar --bzip2 -xf boost_1_59_0.tar.bz2
cd boost_1_59_0
mkdir opt
./bootstrap.sh --prefix=/data/tusers/zhongrenhu/Software/boost_1_59_0/opt
./b2 install

# 创建 LACHESIS 所需目录
mkdir -p opt/stage
ln -s /data/tusers/zhongrenhu/Software/boost_1_59_0/opt/lib /data/tusers/zhongrenhu/Software/boost_1_59_0/opt/stage/lib
```

### 下载LACHESIS 并解压
```
curl -o LACHESIS.zip https://codeload.github.com/shendurelab/LACHESIS/legacy.zip/master
unzip LACHESIS.zip
mv shendurelab-LACHESIS-2e27abb LACHESIS
cd LACHESIS
```

### 修改文件
```
# LACHESIS/src/include/gtools/ 目录下的 SAMStepper.h & SAMStepper.cc
# 对两个文件进行相同的修改

#include <bam/sam.h> ----> #include </data/tusers/zhongrenhu/Software/samtools-0.1.18/sam.h>
```

### 声明环境变量 LD_LIBRARY_PATH
```
将 /data/tusers.ds/zhongrenhu/Software/boost_1_59_0/opt/lib/ 添加到环境变量 LD_LIBRARY_PATH 中，并且在 ~/.profile 文件中进行声明
```

### 声明环境变量 LACHESIS_BOOST_DIR 以及 LACHESIS_SAMTOOLS_DIR
```
# LACHESIS_BOOST_DIR 是 boost 的安装路径
export LACHESIS_BOOST_DIR=/data/tusers/zhongrenhu/Software/boost_1_59_0/opt/

# LACHESIS_SAMTOOLS_DIR 也是 boost 的安装路径(它直接在原本的目录下编译了)
export LACHESIS_SAMTOOLS_DIR=/data/tusers/zhongrenhu/Software/samtools-0.1.18
```
  
### 运行 configure 生成 Makefile
```
./configure --with-samtools=/data/tusers/zhongrenhu/Software/samtools-0.1.18 --with-boost=/data/tusers.ds/zhongrenhu/Software/boost_1_59_0/opt
```

### 修改 LACHESIS/src/include/gtools/Makefile
```
INCLUDES= ----> INCLUDES= -I$(LACHESIS_SAMTOOLS_DIR) -I$(LACHESIS_BOOST_DIR)/include
```

### 修改 LACHESIS/src/include/markov/Makefile
```
INCLUDES= ----> INCLUDES= -I$(LACHESIS_SAMTOOLS_DIR) -I$(LACHESIS_BOOST_DIR)/include

BOOST_LIBS= -lboost_system -lboost_filesystem -lboost_regex ----> BOOST_LIBS = -L$(LACHESIS_BOOST_DIR)/stage/lib -lboost_system -lboost_filesystem -lboost_regex
```

### 运行 make 进行编译
```
make
```

### 使用测试数据尝试运行
```

```