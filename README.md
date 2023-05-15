# STvis
install via:<br>
devtools::install_github("EddieLv/STvis/STvis")

###########################################Tutorial###########################################<br>
library(STvis)<br>
shiny_st()<br>

先看顶部的Instructions熟悉主要内容<br>
软件主要功能有:<br>
-将数值变量、字符变量点亮至切片上<br>
-对字符变量进行筛选并可视化<br>
-自由调整spots(水平、垂直平移；水平、垂直翻转；水平、垂直缩放；360度旋转)直至与切片完全对齐<br>
-基于导入的基因表达矩阵(行与metadata一致, 列是感兴趣的基因)进行可视化<br>
-使用图片右上方的拉索工具进行自由切割(比如可以人工划分出大脑的不同皮层区域)<br>
-下载修改后的metadata文件, 进行后续高级分析(比如基于划分出的大脑的不同皮层区域进行差异分析)

可视化字符变量[small.type]<br>
![image](https://github.com/EddieLv/STvis/assets/61786787/2cfe653c-db25-42a3-ae86-9637cca67a64)
使用图片右上方的[lasso selection]圈出某一个区域<br>
![image](https://github.com/EddieLv/STvis/assets/61786787/02c79294-2b57-4238-bc3a-d5de1bb99751)
在[Set label for selected spots]输入设置的标签, 并点击[Confirm]确认修改<br>
![image](https://github.com/EddieLv/STvis/assets/61786787/b97b4412-7e5c-42ee-94fb-7e0c38e35246)

