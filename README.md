# grtm_cuda
GPU加速的层状介质格林函数计算。

## 快速开始
1. 安装CUDA和cuDNN。
2. 安装依赖包，其实没有啥优化的地方。
3. 安装cmake环境。

编译计算：

```bash
cd build 
cmake ../
make 
./test 
cd ..
python plotfig.py 
python plotfig2.py 
```

