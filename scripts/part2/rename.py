# Pythono3 code to rename multiple
# files in a directory or folder
 
# importing os module
import os
 
# Function to rename multiple files
def main():
   
    folder = "../../data/model1/rename/"
    for filename in os.listdir(folder):
        src = f"{filename}"  # foldername/filename, if .py file is outside folder
        print(src)

        dst = src.split('0.')
        dst = dst[0] + dst[1]

        dst = dst.split('-T0')
        dst = dst[0] + "0-T0" + dst[1]

        if len(dst.split('.0')) == 3:
            print(dst)
            dst = dst.split('.0')
            print(dst)
            dst = dst[0] + dst[1] + dst[2]
            print(dst)

        if len(dst.split('l-')) == 2:
            _, dst = dst.split('l-')
            dst = f"{dst}-optimal"

        if len(dst.split('t-')) == 2:
            _, dst = dst.split('t-')
            dst = f"{dst}-min"
        
        if len(dst.split('x-')) == 2:
            _, dst = dst.split('x-')
            dst = f"{dst}-max"

        dst = f"{folder}/{dst}"

         
        # rename() function will
        # rename all the files
        os.rename(folder+src, dst)


    # for filename in os.listdir(folder):
    #     src = f"{filename}"  # foldername/filename, if .py file is outside folder
    #     os.remove(folder+src)
 
# Driver Code
if __name__ == '__main__':
     
    # Calling main() function
    main()