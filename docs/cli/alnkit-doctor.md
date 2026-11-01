# `alnkit doctor`

打印版本与可用性；必需项缺失时以退出码 1 结束。

## 用法

```bash
alnkit doctor [全局选项]
```

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```bash
alnkit doctor
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
