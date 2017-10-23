
import xlrd
#encoding:utf-8
def open_excel(file='g1vsg2sigregion.xls'):
    try:
        data=xlrd.open_workbook(file)
        return data
    except Exception as e:
        print(str(e))
data=open_excel()
sheet_list=data.sheet_names()   # 获取xls文件中所有sheet的名称
#print(sheet_list)
w=open("result.txt",'w')
w.write("Region\tsum_num\thypomethylated_num\thypermethylated_num\n")
w.close()
def get_summary(sheetn):
    table=data.sheet_by_name(sheetn)  # 通过名称获取一个工作表
    beta_dif=table.col_values(3)[1:]  # 获取整列的值(数组)
    sum=len(beta_dif)
    tf=[i>0 for i in beta_dif]
    w=open("result.txt",'a')
    w.write("%s\t%s\t%s\t%s\n"%(sheetn,sum,tf.count(False),tf.count(True)))
    w.close()
	
for i in range(len(sheet_list)):
    get_summary(sheet_list[i])