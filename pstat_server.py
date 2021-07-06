import socket
import sys


def start_tcp_server(ip, port):
    # create socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_address = (ip, port)

    # bind port
    print("pstat-server: starting listen on ip %s, port %s" % server_address)
    sock.bind(server_address)

    try:
        sock.listen(64)
    except socket.error:
        print("Fail to listen on port %s" % e)
        sys.exit(1)

    connection_id = 0
    is_finish = False
    while True:
        client, addr = sock.accept()
        connection_id += 1
        print("Connection id %d" % connection_id)
        msg = 'welcome to tcp server' + "\r\n"
        receive_count = 0
        receive_count += 1
        while True:
            msg = client.recv(16384)
            msg_de = msg.decode('utf-8')
            print("recv len is : [%d]" % len(msg_de))
            print(msg_de)
            if msg_de == 'disconnect':
                break
            if msg_de == 'server_finish':
                is_finish = True
                break
            msg = ("hello, client, i got your msg %d times, now i will send back to you " % receive_count)
            client.send(msg.encode('utf-8'))
            receive_count += 1
            print("send len is : [%d]" % len(msg))
        if is_finish:
            break
    client.close()
    sock.close()
    print("pstat-server: finished")


if __name__=='__main__':
    start_tcp_server('192.168.1.23', 6000)
