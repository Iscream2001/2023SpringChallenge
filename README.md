# ICPC 2023 Online Spring Challenge powered by Huawei 记录


## 整体思路

我的思路：

1. 在每个用户内先选出一个最应该被丢弃的。

2. 然后再比较不同用户谁的最应该被替换。


## 同一用户的内存丢弃

### 丢弃规则（比较方式）

这里我使用3（4）种丢弃的比较方式：

1. 丢弃掉最久未被使用的

2. 丢弃掉被使用次数最多的

3. 丢弃掉被使用次数最少的

4. 丢弃掉最近被使用的 （但是因为原数据集上似乎没有适合这种丢弃方法的数据，并且我的代码在原数据集上跑得很慢，我最后把这种方法去掉了，然后system test爆零了一个点）。

### 比较方式的选取

以上几种方法在不同类型的数据（或者同一数据的不同时期）上效果差异很大。

不同的情况需要选取不同的比较方式。

我是这样选择的：

1. 对于每种比较方式，新增一个与真实算法同步的一个List，然后这个List的每个用户体积与真实算法的体积分配情况相同，然后使用相应的比较方法，记录命中失败的次数。

2. 比较3（4）个同步List的命中失败次数，把失败次数最少的策略应用到真是算法的下一时刻。

## 不同用户间的内存丢弃

关于这个我没有想到什么特别好的想法。

我最后的做法：

1. 考虑与用户间的丢弃相关的因素：用户的优先级，本次丢弃导致的代价函数上升的幅度，内存本身的因素（上次使用时间，使用频率等）。

2. 考虑将三种因素的值乘起来，作为丢弃这个内存会产生的代价。

3. 由于时间跨度很长，而且用户被丢弃一次，其占用的总内存数就会下降，猜测丢弃就必然导致其代价会被计入最终的总代价。所以总是选取代价小的内存来丢弃。

4. 代价的3个因素，都可以通过调整指数的方式进行调参（我代码里那些奇奇怪怪的参数就是这么来的）。


## 比赛结果

最终跑出来的结果跟我预想的还是差不多的，几乎所有数据点都在400-500分的亚子。

只有某个点爆零了（估计如果当时不删第4个比较方式可能就不会爆零了，这个点如果那个300或者400就能前10惹QAQ）

最后成绩是39965，rank26 。

