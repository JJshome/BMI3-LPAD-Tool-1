import os

def read_file(file_path):
    """读取文件内容"""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"文件 {file_path} 不存在")
    with open(file_path, 'r') as file:
        return file.read()

def write_file(file_path, content):
    """将内容写入文件"""
    with open(file_path, 'w') as file:
        file.write(content)
